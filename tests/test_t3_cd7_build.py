"""Row R36's file cD7 in experiments/t3_decomp/build_two_knockout.py: C2 on the barred s2 (`--form c2-s2`) rebuilt by
its own command with one change, its composition block replaced by the decision file's
submissions.allowed.s1.composition (`--composition-from allowed-s1`, only with `--form c2-s2` and `--captain-ruling
R36`, only while row R36's first cell names it), and its two target-free gates: g1 (the composition routine finds the
five scenario states at 0.3835 of their wild-type share before renormalising) and g2 (the file's gate bar, as row N4's
diagnostic computes it, at most C2's plus 0.01).

Offline, synthetic and fast: the throwaway world of tests/test_t3_two_knockout_build.py (a small synthetic WT file,
synthetic `transfer_mab21l2_all` and `gata4_pooled` tables, a decision file with the frozen hashes), a decision file
that also carries the allowed s1 composition block, and small ledger files written under tmp_path. Nothing here reads
the repository's own TRIES.md, builds from released data or pins bytes built on one machine.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
T3 = ROOT / "experiments" / "t3_decomp"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


th = sys.modules.get("_t3_two_knockout_helpers") or _load("_t3_two_knockout_helpers",
                                                          ROOT / "tests" / "test_t3_two_knockout_build.py")
bk, bs, tx = th.bk, th.bs, th.tx
t3d = bs.t3d
world, synthetic_pins, own_repo = th.world, th.synthetic_pins, th.own_repo  # fixtures (the last two autouse)
SEED, N, G, GENES = th.SEED, th.N, th.G, th.GENES
REAL_DECISION_PIN = bk.R36_DECISION_SHA256  # read before any test patches it
RULES = tx.RULES
STATES = list(RULES["presets"]["gata4"]["scenario"]["depleted_states"])
L2 = 1.3826521560990648
# the block as the frozen rules write a central composition at this magnitude (decide_submissions' form)
BLOCK = json.loads(json.dumps(bs.composition_spec("central", L2, RULES, "Gata4 replicates (D7)")))
OPT = "`--composition-from allowed-s1`"
ROW = ("| R36 (the captain's named exception of a day, \"yes, one file\": a composition changes for one file; "
       "ruling, no run, not a try; the T3 two-knockout builder is to take `--composition-from allowed-s1` with "
       "`--form c2-s2` only with `--captain-ruling R36`) | T3 | x |\n")
OTHER_FORMS = [f for f in ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2",
                           "c2-s1", "c2-s1x1.5", "c2-s1x2.0")]


def write_decision(path: Path, block=BLOCK, **over) -> Path:
    """A synthetic decision file: the frozen hashes, D7m, and submissions.allowed.s1.composition = block."""
    sub = {"allowed": {"s1": {"composition": block}}} if block is not None else {}
    return tx._decision(path, **({"submissions": sub} | over))


@pytest.fixture(scope="module")
def cw(world):
    dec = write_decision(world.d / "gata4_decisions_cd7.json")
    return SimpleNamespace(dec=dec, sha=t3d.sha256(dec))


@pytest.fixture(autouse=True)
def decision_pin(cw, monkeypatch):
    """The builder pins the committed decision file's sha256; here the synthetic one's (the pin itself is tested)."""
    monkeypatch.setattr(bk, "R36_DECISION_SHA256", cw.sha)


def ledger(tmp_path: Path, rows: str = ROW, name: str = "TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text(f"# a throwaway ledger\n| R22 (another row) | T3 | y |\n{rows}| R37 (a later row) | T1 | z |\n")
    return p


def build(world, cw, out: Path, form: str = "c2-s2", **kw) -> dict:
    return bk.build_form(kw.pop("dec", cw.dec), form, kw.pop("wt", world.wt), world.board, N, out,
                         run_dir=kw.pop("run_dir", world.run), seed=SEED, threads=2, **kw)


def build_cd7(world, cw, tmp_path: Path, out: Path, **kw) -> dict:
    return build(world, cw, out, composition_from="allowed-s1", captain_ruling="R36", tries=ledger(tmp_path), **kw)


def nothing_written(out: Path) -> bool:
    return not out.exists() and not out.with_suffix(".json").exists() and not Path(f"{out}.stamp.json").exists()


def frozen_functions_restored() -> bool:
    return all(getattr(bs, n).__module__ == "build_submission" and getattr(bs, n).__name__ == n
               for n in ("resolve_spec", "composition", "load_factor"))


class Spy:
    """Records the factor apply_gene_factor receives."""

    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


def wt_without(tmp_path: Path, monkeypatch, state: str) -> SimpleNamespace:
    """A synthetic WT stage without one scenario state, and a run folder whose record names it (same tables)."""
    with monkeypatch.context() as m:
        m.setattr(tx, "CELLS", [c for c in tx.CELLS if state not in (c[0], c[1])])
        wt = tmp_path / f"wt_no_{state.replace('/', '_')}.h5ad"
        tx._write_wt(wt, np.random.default_rng(1))
    run = tmp_path / f"run_no_{state.replace('/', '_')}"
    th.write_run(run, wt, {bk.TABLE_M: th.PHI_M, bk.TABLE_G: th.PHI_G})
    return SimpleNamespace(wt=wt, run=run)


# ------------------------------------------------------------------------------------------------------ constants


def test_the_line_s_constants():
    assert (bk.R36_ROW, bk.R36_FORM, bk.R36_BASE, bk.R36_OPTION, bk.R36_FILE) == (
        "R36", "c2-s2", "s2", OPT, "T3_c2_s2_cd7.h5ad")
    assert bk.R36_NAMES == (OPT, "`--form c2-s2`", "`--captain-ruling R36`")
    assert bk.R36_OPENER == "the captain's named exception"
    assert bk.COMPOSITION_FROM == {"allowed-s1": ("allowed", "s1")}
    assert bk.R36_LOG2_MAGNITUDE == L2 and bk.R36_RATIO == 2.0 ** -L2 and round(bk.R36_RATIO, 4) == 0.3835
    assert bk.R36_N_STATES == 5 == len(STATES) and bk.R36_G2_MARGIN == 0.01 and bk.R36_RATIO_RTOL == 1e-9
    assert (bk.EXIT_REFUSED, bk.EXIT_GATE_FAILED, bk.EXIT_GATE_CRASHED, bk.EXIT_FACTOR_MISMATCH) == (2, 1, 3, 5)
    assert bk.N4_SCRIPT == ROOT / "experiments" / "research_2026_11" / "diagnostics.py" and bk.N4_SCRIPT.is_file()
    assert bk.TRIES_FILE == ROOT / "experiments" / "harness" / "TRIES.md"
    # C2's recorded sha256 is the build record's
    words = (ROOT / "experiments" / "builds" / "results" / "t3two_fb01487.SHA256SUMS").read_text().split()
    assert dict(zip(words[1::2], words[0::2], strict=True))["c2-s2.h5ad"] == bk.C2_S2_SHA256
    # the line adds no form and changes no kind, base or row R22 constant
    assert bk.FORMS == ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2",
                        "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2")
    assert bk.KINDS == {"identity": 0.0, "c1": 1.0, "c2": 1.0} and bk.BASES["s2"][:2] == ("s2", 1.0)
    assert (bk.G2_ROW, bk.G2_FORM, bk.G2_OPTION, bk.G2_FILE) == ("R22", "c1-s2", "`--gata4-exponent 2.0`",
                                                                 "T3_c1_s2_g2.h5ad")
    # build_outside_form's helpers are reused there, never defined here (row R22's own test holds this too)
    assert not hasattr(bk, "captain_order_names") and not hasattr(bk, "wt_den")


def test_the_pin_is_the_committed_decision_file_s_and_its_block_is_the_line_s(monkeypatch):
    """The pin is the committed decision file's sha256, and that file's submissions.allowed.s1.composition passes the
    builder's check: the frozen rules' central composition at the pinned magnitude, as composition_spec writes it."""
    path = T3 / "results" / "gata4" / "gata4_decisions.json"
    assert REAL_DECISION_PIN == hashlib.sha256(path.read_bytes()).hexdigest()
    monkeypatch.setattr(bk, "R36_DECISION_SHA256", REAL_DECISION_PIN)
    dec = json.loads(path.read_text())
    block, rec = bk.composition_block(dec, RULES, path)
    assert block == dec["submissions"]["allowed"]["s1"]["composition"] == BLOCK
    assert block["depleted_states"] == STATES and block["log2_magnitude"] == L2
    assert block["depletion_factor"] == bk.R36_RATIO
    assert rec["decision_sha256"] == rec["pinned_sha256"] == REAL_DECISION_PIN
    assert rec["source"] == "decision file submissions.allowed.s1.composition"


# ------------------------------------------------------------------------------------------------- the ruling row


@pytest.mark.parametrize("rows,ok", [
    (ROW, True),
    (ROW + ROW, False),  # two rows R36
    ("", False),  # no row R36
    (ROW.replace("| R36 (the captain's named exception", "| R36 (the captain's order"), False),
    (ROW.replace("| R36 (the captain's named exception", "| R36 (a named exception"), False),
    (ROW.replace("ruling, no run", "a run"), False),
    (ROW.replace(OPT, "`--composition-from allowed-s2`"), False),
    (ROW.replace(OPT, "--composition-from allowed-s1"), False),  # not in backticks
    (ROW.replace("`--form c2-s2`", "`--form c1-s2`"), False),
    (ROW.replace("`--captain-ruling R36`", "`--captain-ruling R22`"), False),
    (ROW.replace(f"take {OPT} with", "take it with").replace("| T3 | x |", f"| T3 | {OPT} |"),
     False),  # the option only outside the first cell
    (ROW.replace("| R36 (", "| R36b ("), False),
])
def test_the_ruling_row_check(tmp_path, rows, ok):
    p = ledger(tmp_path, rows)
    assert bk.ruling_row_names("R36", bk.R36_OPENER, bk.R36_NAMES, p) is ok
    if ok:
        rec = bk.check_composition_from("c2-s2", "allowed-s1", "R36", None, p)
        assert rec == {"row": "R36", "option": OPT, "composition_from": "allowed-s1",
                       "ledger": "experiments/harness/TRIES.md", "ledger_sha256": t3d.sha256(p)}
    else:
        with pytest.raises(bk.Refused, match="no single row R36"):
            bk.check_composition_from("c2-s2", "allowed-s1", "R36", None, p)


def test_a_missing_ledger_is_refused(tmp_path):
    with pytest.raises(bk.Refused, match="no single row R36"):
        bk.check_composition_from("c2-s2", "allowed-s1", "R36", None, tmp_path / "absent.md")


def test_without_the_option_or_the_ruling_r36_nothing_is_read(monkeypatch):
    monkeypatch.setattr(bk, "ruling_row_names", lambda *a, **k: pytest.fail("the ledger was read"))
    for form in [*bk.FORMS, "c1-s1x1.5", "c3-s1", "", None]:
        for kw in ({}, {"captain_ruling": "R22"}, {"captain_ruling": "R22", "gata4_exponent": 2.0},
                   {"gata4_exponent": 2.0}, {"captain_ruling": "R30"}, {"captain_ruling": "R36 "},
                   {"captain_ruling": "r36"}):
            assert bk.check_composition_from(form, **kw) is None


# ------------------------------------------------------------------------------------------------------ refusals


@pytest.mark.parametrize("value", ["allowed-s2", "refused-s1", "barred-s1", "allowed_s1", "Allowed-s1",
                                   "allowed-s1 ", " allowed-s1", "", "s1", "allowed", 1, True, ["allowed-s1"]])
def test_any_other_value_is_refused(tmp_path, value):
    with pytest.raises(bk.Refused, match="only allowed-s1"):
        bk.check_composition_from("c2-s2", value, "R36", None, ledger(tmp_path))


@pytest.mark.parametrize("form", [*OTHER_FORMS, "c2-s2 ", "C2-s2", None])
def test_any_other_form_is_refused(tmp_path, form):
    with pytest.raises(bk.Refused, match="on --form c2-s2 only"):
        bk.check_composition_from(form, "allowed-s1", "R36", None, ledger(tmp_path))


@pytest.mark.parametrize("ruling", [None, "R22", "R30", "r36", "R36 ", " R36", "R36b", "", "K23"])
def test_the_option_without_the_ruling_r36_is_refused(tmp_path, ruling):
    with pytest.raises(bk.Refused, match="only with --captain-ruling R36"):
        bk.check_composition_from("c2-s2", "allowed-s1", ruling, None, ledger(tmp_path))


def test_the_ruling_without_the_option_and_the_option_with_another_line_s(tmp_path):
    for form in ("c2-s2", "c1-s2", "identity-s2"):
        for g in (None, 2.0):
            with pytest.raises(bk.Refused, match="taken only with --composition-from allowed-s1"):
                bk.check_composition_from(form, None, "R36", g, ledger(tmp_path))
    for g in (2.0, 1.0, "2.0"):
        with pytest.raises(bk.Refused, match="takes no --gata4-exponent"):
            bk.check_composition_from("c2-s2", "allowed-s1", "R36", g, ledger(tmp_path))


@pytest.mark.parametrize("kw", [
    {"composition_from": "allowed-s2", "captain_ruling": "R36"},
    {"composition_from": "allowed-s1"},
    {"captain_ruling": "R36"},
    {"composition_from": "allowed-s1", "captain_ruling": "R22"},
    {"composition_from": "allowed-s1", "captain_ruling": "R36", "form": "c1-s2"},
    {"composition_from": "allowed-s1", "captain_ruling": "R36", "form": "identity-s2"},
    {"composition_from": "allowed-s1", "captain_ruling": "R36", "gata4_exponent": 2.0},
    {"composition_from": "allowed-s1", "captain_ruling": "R36", "rows": ""},
    {"composition_from": "allowed-s1", "captain_ruling": "R36", "rows": ROW.replace(OPT, "`--composition-from`")},
    {"captain_ruling": "R36", "gata4_exponent": 2.0, "form": "c1-s2"},
])
def test_the_builder_refuses_and_builds_nothing(world, cw, tmp_path, kw, monkeypatch):
    for name in ("composition_block", "r36_composition", "run_gata4_figures"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on a refusal"))
    kw = dict(kw)
    form = kw.pop("form", "c2-s2")
    t = ledger(tmp_path, kw.pop("rows", ROW))
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused):
        build(world, cw, out, form, tries=t, **kw)
    assert nothing_written(out) and frozen_functions_restored()


@pytest.mark.parametrize("block,says", [
    ({**BLOCK, "log2_magnitude": 0.7}, "not the central composition"),
    ({**BLOCK, "log2_magnitude": float(np.nextafter(L2, 0.0))}, "not the central composition"),  # one ulp below
    ({**BLOCK, "log2_magnitude": L2 / 2}, "not the central composition"),  # the mild magnitude of the same m
    ({**BLOCK, "log2_magnitude": "1.3826521560990648"}, "not the central composition"),
    ({**BLOCK, "kind": "mild"}, "not the central composition"),
    ({**BLOCK, "depleted_states": STATES[::-1]}, "not the central composition"),
    ({**BLOCK, "depleted_states": STATES[:4]}, "not the central composition"),
    ({**BLOCK, "spared_states": BLOCK["spared_states"][:2]}, "not the central composition"),
    ({**BLOCK, "depletion_factor": 0.3835}, "not the central composition"),
    ({k: v for k, v in BLOCK.items() if k != "depletion_factor"}, "not the central composition"),
    ("central", "not the central composition"),
    (None, "has no submissions.allowed.s1.composition"),
])
def test_a_block_other_than_the_line_s_is_refused(world, cw, tmp_path, monkeypatch, block, says):
    dec = write_decision(tmp_path / "dec.json", block)
    monkeypatch.setattr(bk, "R36_DECISION_SHA256", t3d.sha256(dec))
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match=says):
        build_cd7(world, cw, tmp_path, out, dec=dec)
    assert nothing_written(out) and frozen_functions_restored()


def test_a_decision_file_other_than_the_pinned_one_is_refused(world, cw, tmp_path, monkeypatch):
    monkeypatch.setattr(bk, "R36_DECISION_SHA256", REAL_DECISION_PIN)  # the synthetic file is not the committed one
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match="not the pinned one"):
        build_cd7(world, cw, tmp_path, out)
    assert nothing_written(out)
    # the same block under another file: refused too (the block is read from the pinned file only)
    other = write_decision(tmp_path / "other.json", seeds=list(range(11)))
    with pytest.raises(bk.Refused, match="not the pinned one"):
        build_cd7(world, cw, tmp_path, out, dec=other)
    assert nothing_written(out)


# ------------------------------------------------------------------------------------------------------- the file


def test_cd7_is_c2_with_its_composition_block_replaced(world, cw, tmp_path, monkeypatch):
    spy = Spy(monkeypatch)
    c2 = build(world, cw, tmp_path / "c2.h5ad")
    out = tmp_path / bk.R36_FILE
    cd7 = build_cd7(world, cw, tmp_path, out)
    # C2's factor, applied once in each
    assert len(spy.seen) == 2 and np.array_equal(spy.seen[0], spy.seen[1])
    # the spec: C2's, its composition block replaced by the decision file's, nothing else
    assert c2["spec"]["composition"]["kind"] == "mild" and cd7["spec"]["composition"] == BLOCK
    assert cd7["spec"] == {**c2["spec"], "composition": BLOCK}
    assert cd7["spec_source"] == {**c2["spec_source"], "composition_source": bk.R36_SOURCE}
    assert cd7["branch"] == c2["branch"] == "barred" and cd7["submission"] == c2["submission"] == "s2"
    # factor, row number, seed, streams, board and coordinates rule: C2's
    assert cd7["inputs"] == c2["inputs"]
    assert cd7["inputs"]["rng_streams"] == ["submission:barred:s2", "submission:barred:s2:assemble"]
    assert cd7["inputs"]["seed"] == SEED and cd7["warnings"] == c2["warnings"]
    con7, con2 = cd7["construction"], c2["construction"]
    assert {k: v for k, v in con7.items() if k not in ("composition", "cells_by_source", "duplicate_positions_moved")
            } == {k: v for k, v in con2.items() if k not in ("composition", "cells_by_source",
                                                             "duplicate_positions_moved")}
    assert con7["factor_exponent"] == 1.0 and con7["depth_ratios"] is None and con7["stage_alphas"] is None
    comp = con7["composition"]
    assert comp["kind"] == "central" and comp["log2_magnitude"] == L2
    assert comp["states_found"] == STATES and comp["states_missing"] == []
    # combined_build: C2's, plus the line's block
    cb7, cb2 = cd7["combined_build"], c2["combined_build"]
    assert {k: v for k, v in cb7.items() if k != "composition_line"} == cb2
    cl = cb7["composition_line"]
    assert cl["composition"] == BLOCK and cl["replaced"] == c2["spec"]["composition"]
    assert cl["captain_ruling"]["row"] == "R36" and cl["captain_ruling"]["option"] == OPT
    assert cl["composition_from"] == "allowed-s1" and cl["file"] == "T3_c2_s2_cd7.h5ad"
    assert cl["decision_file"]["decision_sha256"] == cw.sha and cl["log2_magnitude"] == L2
    assert list(cl["ratios_before_renormalising"]) == STATES
    assert all(r == pytest.approx(bk.R36_RATIO, rel=1e-12, abs=0) for r in cl["ratios_before_renormalising"].values())
    assert cl["c2_reference_sha256"] == bk.C2_S2_SHA256 and cl["line"] == bk.R36_LINE
    assert not any("cell" in k.lower() for k in cl)
    on_disk = json.loads(out.with_suffix(".json").read_text())
    assert on_disk["combined_build"]["composition_line"]["composition"] == BLOCK and on_disk["spec"] == cd7["spec"]
    # a fresh draw (renormalising moves every quota), format-checked, stamped clean
    assert cd7["output"]["check"]["ok"] and (tmp_path / "c2.h5ad").read_bytes() != out.read_bytes()
    stamp = json.loads(Path(f"{out}.stamp.json").read_text())
    assert stamp["clean"] is True and stamp["sha256"] == cd7["output"]["sha256"] == t3d.sha256(out)
    assert frozen_functions_restored()


def test_the_shares_are_the_scenario_at_0_3835_by_hand(world, cw, tmp_path):
    side = build_cd7(world, cw, tmp_path, tmp_path / "f.h5ad")
    ref = bs.load_reference(world.wt, "gata4", RULES, None, 2)
    p = np.asarray(ref.V.p, dtype=np.float64)
    raw = p.copy()
    idx = ref.V.index(STATES)
    assert len(idx) == 5
    raw[idx] *= 2.0 ** -L2
    want = raw / raw.sum()
    shares = side["construction"]["composition"]["shares"]
    got = np.array([shares[t] for t in ref.V.labels])
    assert np.allclose(got, want, rtol=1e-12, atol=0)
    rest = np.setdiff1d(np.arange(p.size), idx)
    assert np.allclose(got[idx] / p[idx] / (got[rest] / p[rest]).mean(), 0.3835131214774091, rtol=1e-12, atol=0)


def test_built_twice_gives_one_sha256(world, cw, tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = build_cd7(world, cw, tmp_path, tmp_path / "a" / bk.R36_FILE)
    b = build_cd7(world, cw, tmp_path, tmp_path / "b" / bk.R36_FILE)
    assert (tmp_path / "a" / bk.R36_FILE).read_bytes() == (tmp_path / "b" / bk.R36_FILE).read_bytes()
    assert a["output"]["sha256"] == b["output"]["sha256"]


def test_the_build_refuses_a_routine_that_does_not_give_the_five_states(world, cw, tmp_path, monkeypatch):
    # a reference without one of the five states: the frozen builder would only warn; this line refuses, unwritten
    other = wt_without(tmp_path, monkeypatch, "pPHM")
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match="does not give the 5 scenario states"):
        build_cd7(world, cw, tmp_path, out, wt=other.wt, run_dir=other.run)
    assert nothing_written(out) and frozen_functions_restored()
    # a routine that moves one state otherwise
    real = bs.composition

    def skewed(V, spec, params):
        f, info = real(V, spec, params)
        f = f.copy()
        f[V.index(["aSHF"])] *= 1.01
        return f / f.sum(), info

    monkeypatch.setattr(bs, "composition", skewed)
    with pytest.raises(bk.Refused, match="does not give the 5 scenario states"):
        build_cd7(world, cw, tmp_path, out)
    assert nothing_written(out)


@pytest.mark.parametrize("form", list(bk.FORMS))
def test_every_existing_form_never_touches_the_line(world, cw, tmp_path, monkeypatch, form):
    """Without the option the line's code is never reached: each form builds the same bytes and sidecar block with
    the line's functions booby-trapped as without, and no sidecar has the line's block."""
    plain = build(world, cw, tmp_path / "plain.h5ad", form)
    for name in ("ruling_row_names", "composition_block", "with_composition", "composition_ratios", "ratios_pass",
                 "r36_composition", "run_g1", "run_g2", "_n4", "n4_bar", "g2_passes"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on {form}"))
    trapped = build(world, cw, tmp_path / "trapped.h5ad", form, gata4_exponent=None, captain_ruling=None,
                    composition_from=None)
    assert (tmp_path / "plain.h5ad").read_bytes() == (tmp_path / "trapped.h5ad").read_bytes()
    assert plain["combined_build"] == {**trapped["combined_build"]}
    assert "composition_line" not in plain["combined_build"] and plain["spec"] == trapped["spec"]
    assert "composition_source" not in plain["spec_source"]


def test_row_r22_s_option_set_never_touches_the_line(world, cw, tmp_path, monkeypatch):
    r22 = tmp_path / "r22.md"
    r22.write_text("| R22 (the captain's order of a day; ruling, no run; the T3 two-knockout builder is to take "
                   "`--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`) | T3 | x |\n")
    plain = build(world, cw, tmp_path / "plain.h5ad", "c1-s2", gata4_exponent=2.0, captain_ruling="R22", tries=r22)
    for name in ("ruling_row_names", "composition_block", "r36_composition", "composition_ratios"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on row R22's line"))
    trapped = build(world, cw, tmp_path / "trapped.h5ad", "c1-s2", gata4_exponent=2.0, captain_ruling="R22",
                    tries=r22, composition_from=None)
    assert (tmp_path / "plain.h5ad").read_bytes() == (tmp_path / "trapped.h5ad").read_bytes()
    assert "composition_line" not in trapped["combined_build"] and "gata4_line" in trapped["combined_build"]
    assert plain["combined_build"]["gata4_line"]["factor_sha256_applied"] == trapped["combined_build"][
        "gata4_line"]["factor_sha256_applied"]


# --------------------------------------------------------------------------------------------------------- G1


def test_composition_ratios_by_hand():
    V = SimpleNamespace(labels=["a", "b", "c", "d", "e"], p=np.array([0.1, 0.2, 0.3, 0.25, 0.15]))
    raw = V.p.copy()
    raw[[0, 1]] *= 0.5
    f = raw / raw.sum()
    info = {"states_found": ["a", "b"], "states_missing": []}
    assert bk.composition_ratios(V, f, info, ["a", "b"]) == pytest.approx([0.5, 0.5], rel=1e-14)
    # a state the routine did not find, or a named state absent from the vocabulary: nan for it
    r = bk.composition_ratios(V, f, {"states_found": ["a"], "states_missing": ["b"]}, ["a", "b"])
    assert r[0] == pytest.approx(0.5) and np.isnan(r[1])
    r = bk.composition_ratios(V, f, info, ["a", "z"])
    assert np.isnan(r[1])
    # the routine changed another state too: no single renormalising factor, nan for every state
    g = raw.copy()
    g[2] *= 1.5
    assert all(np.isnan(bk.composition_ratios(V, g / g.sum(), info, ["a", "b"])))
    # no state left alone
    assert all(np.isnan(bk.composition_ratios(V, f, info, ["a", "b", "c", "d", "e"])))


def test_ratios_pass():
    r = bk.R36_RATIO
    assert bk.ratios_pass([r] * 5) and bk.ratios_pass([r * (1 + 5e-10)] * 5)
    assert not bk.ratios_pass([r] * 4) and not bk.ratios_pass([r] * 6)
    assert not bk.ratios_pass([r, r, r, r, r * (1 + 1e-6)]) and not bk.ratios_pass([r, r, r, r, float("nan")])
    assert not bk.ratios_pass([0.3835] * 5)  # the rounded value is not the routine's
    assert not bk.ratios_pass([2.0 ** (-L2 / 2)] * 5)  # the mild magnitude


def test_g1_finds_the_five_states_at_0_3835(world, cw):
    ratios = bk.run_g1(cw.dec, world.wt, "allowed-s1", threads=2)
    assert len(ratios) == 5 and all(x == pytest.approx(bk.R36_RATIO, rel=1e-12, abs=0) for x in ratios)
    assert bk.ratios_pass(ratios) and frozen_functions_restored()


def _g1(args):
    from typer.testing import CliRunner

    return CliRunner().invoke(bk.app, ["g1", *map(str, args)])


def test_g1_cli_prints_the_five_ratios_and_nothing_else(world, cw):
    res = _g1([cw.dec, "--wt", world.wt, "--composition-from", "allowed-s1", "--threads", "2"])
    assert res.exit_code == 0, res.output
    assert res.stdout.splitlines() == ["0.3835"] * 5


def test_g1_fails_without_a_state_or_with_a_routine_that_moves_one_otherwise(world, cw, tmp_path, monkeypatch):
    other = wt_without(tmp_path, monkeypatch, "pPHM")
    res = _g1([cw.dec, "--wt", other.wt, "--composition-from", "allowed-s1", "--threads", "2"])
    assert res.exit_code == 1, res.output
    assert res.stdout.splitlines() == ["0.3835"] * 4 + ["nan"]  # pPHM is the block's fifth state
    real = bs.composition

    def skewed(V, spec, params):
        f, info = real(V, spec, params)
        f = f.copy()
        f[V.index(["aSHF"])] *= 1.01
        return f / f.sum(), info

    monkeypatch.setattr(bs, "composition", skewed)
    res = _g1([cw.dec, "--wt", world.wt, "--composition-from", "allowed-s1", "--threads", "2"])
    assert res.exit_code == 1 and len(res.stdout.splitlines()) == 5


def test_g1_refusals_and_a_failure_print_no_figure(world, cw, monkeypatch):
    for value in ("allowed-s2", "refused-s1", ""):
        res = _g1([cw.dec, "--wt", world.wt, "--composition-from", value, "--threads", "2"])
        assert res.exit_code == 2 and "refused" in res.output
    monkeypatch.setattr(bk, "R36_DECISION_SHA256", REAL_DECISION_PIN)
    res = _g1([cw.dec, "--wt", world.wt, "--composition-from", "allowed-s1", "--threads", "2"])
    assert res.exit_code == 2 and "not the pinned one" in res.output
    monkeypatch.setattr(bk, "run_g1", lambda *a, **k: (_ for _ in ()).throw(ValueError("shape (2500, 40)")))
    res = _g1([cw.dec, "--wt", world.wt, "--composition-from", "allowed-s1"])
    assert res.exit_code == 3 and res.stdout.strip() == "failed: ValueError" and "2500" not in res.output


# --------------------------------------------------------------------------------------------------------- G2


def test_n4_bar_is_row_n4_s_gate_bar():
    n4 = bk._n4()
    assert Path(n4.__file__).resolve() == bk.N4_SCRIPT.resolve() and n4 is bk._n4()
    rng = np.random.default_rng(3)
    dp = rng.normal(0.01, 0.1, 500)
    assert bk.n4_bar(dp) == 0.01 + 0.99 * n4.cos2_ones(dp)
    assert bk.n4_bar(dp) == pytest.approx(0.01 + 0.99 * dp.sum() ** 2 / (dp.size * float(dp @ dp)), rel=1e-12)
    assert bk.n4_bar(np.full(7, -0.2)) == pytest.approx(1.0) and bk.n4_bar(np.array([1.0, -1.0])) == 0.01


def test_g2_passes():
    assert bk.g2_passes(0.10, 0.10) and bk.g2_passes(0.10, 0.05) and bk.g2_passes(0.10, 0.1099)
    assert not bk.g2_passes(0.10, 0.1101) and not bk.g2_passes(0.02, 0.5)
    assert not bk.g2_passes(float("nan"), 0.1) and not bk.g2_passes(0.1, float("nan"))


def test_g2_bars_are_the_two_built_files_bars(world, cw, tmp_path):
    """G2 draws each file's cells from its own build inputs, in memory; its bars are the bars of the files the
    builder writes (pseudobulk over the board's panel minus WT's as released, N4's functions), to float32 summation."""
    import anndata as ad

    from vec.prep import load_stage

    res = bk.run_g2(cw.dec, world.wt, world.run, world.board, N, SEED, "allowed-s1", threads=2)
    assert frozen_functions_restored()
    assert set(res) == {"c2", "cd7", "margin", "gene_set", "passes"} and res["margin"] == 0.01
    assert res["passes"] == bk.g2_passes(res["c2"], res["cd7"]) and res["gene_set"] == "all panel genes"
    n4 = bk._n4()
    w = load_stage(world.wt)
    pos = {g: i for i, g in enumerate(w.var_names)}
    pb_wt = n4.pseudobulk(w.X.tocsr()[:, [pos[g] for g in world.board.genes]])
    build(world, cw, tmp_path / "c2.h5ad")
    build_cd7(world, cw, tmp_path, tmp_path / "cd7.h5ad")
    for name in ("c2", "cd7"):
        a = ad.read_h5ad(tmp_path / f"{name}.h5ad")
        assert list(a.var_names) == list(world.board.genes)
        dp = n4.pseudobulk(a.X) - pb_wt
        assert res[name] == pytest.approx(bk.n4_bar(dp), abs=1e-6)
    assert res["c2"] != res["cd7"]


def _g2(args):
    from typer.testing import CliRunner

    return CliRunner().invoke(bk.app, ["g2", *map(str, args)])


def _g2_args(world, cw, tmp_path, *extra):
    return [cw.dec, "--wt", world.wt, "--board", "T3:synthetic", "--boards-dir", th.cli_boards(tmp_path),
            "--run-dir", world.run, "--composition-from", "allowed-s1", "--threads", "2", *extra]


def test_g2_cli_prints_the_two_bars_and_nothing_else(world, cw, tmp_path):
    res = _g2(_g2_args(world, cw, tmp_path))
    want = bk.run_g2(cw.dec, world.wt, world.run, world.board, N, SEED, "allowed-s1", threads=2)
    assert res.exit_code == (0 if want["passes"] else 1), res.output
    assert res.stdout.splitlines() == [f"{want['c2']:.4f}", f"{want['cd7']:.4f}"]


@pytest.mark.parametrize("c2,cd7,code", [(0.10, 0.105, 0), (0.10, 0.11, 0), (0.10, 0.2, 1), (0.30, 0.01, 0)])
def test_g2_cli_exits_by_the_verdict(world, cw, tmp_path, monkeypatch, c2, cd7, code):
    seen = []

    def fake(decision, wt, run_dir, board, n_cells, seed, composition_from, threads):
        seen.append((n_cells, seed, composition_from))
        return {"c2": c2, "cd7": cd7, "margin": 0.01, "gene_set": "all panel genes", "passes": bk.g2_passes(c2, cd7)}

    monkeypatch.setattr(bk, "run_g2", fake)
    res = _g2(_g2_args(world, cw, tmp_path))
    assert res.exit_code == code and res.stdout.splitlines() == [f"{c2:.4f}", f"{cd7:.4f}"]
    assert seen == [(N, SEED, "allowed-s1")]  # the board's max_cells, read in code


def test_g2_refusals_and_a_failure_print_no_figure(world, cw, tmp_path, monkeypatch):
    boards = th.cli_boards(tmp_path)
    base = [cw.dec, "--wt", world.wt, "--board", "T3:synthetic", "--boards-dir", boards, "--run-dir", world.run,
            "--threads", "2"]
    for extra in (["--composition-from", "allowed-s2"], ["--composition-from", "allowed-s1", "--seed", "1000002"],
                  ["--composition-from", "allowed-s1", "--n-cells", str(N)], ["--composition-from", "", ]):
        res = _g2([*base, *extra])
        assert res.exit_code == 2 and "refused" in res.output, (extra, res.output)
    monkeypatch.setattr(bk, "R36_DECISION_SHA256", REAL_DECISION_PIN)
    res = _g2([*base, "--composition-from", "allowed-s1"])
    assert res.exit_code == 2 and "not the pinned one" in res.output
    monkeypatch.setattr(bk, "run_g2", lambda *a, **k: (_ for _ in ()).throw(ValueError(f"{N} cells")))
    res = _g2([*base, "--composition-from", "allowed-s1"])
    assert res.exit_code == 3 and res.stdout.strip() == "failed: ValueError" and str(N) not in res.output


# --------------------------------------------------------------------------------------------------------- CLI


def test_the_build_cli(world, cw, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    from vec.io import check_submission

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    base = ["build", str(cw.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(th.cli_boards(tmp_path)), "--run-dir", str(world.run), "--seed", str(SEED), "--threads", "2"]
    runner = CliRunner()
    for extra in (["--form", "c2-s2", "--composition-from", "allowed-s2", "--captain-ruling", "R36"],
                  ["--form", "c2-s2", "--composition-from", "allowed-s1"],
                  ["--form", "c2-s2", "--captain-ruling", "R36"],
                  ["--form", "c1-s2", "--composition-from", "allowed-s1", "--captain-ruling", "R36"],
                  ["--form", "c2-s2", "--composition-from", "allowed-s1", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--composition-from", "allowed-s1", "--captain-ruling", "R36",
                   "--gata4-exponent", "2.0"],
                  ["--form", "c1-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R36"],
                  ["--form", "c2-s2", "--composition-from", "allowed-s1", "--captain-ruling", "R36", "--seed",
                   "1000002"]):
        out = tmp_path / "r.h5ad"
        res = runner.invoke(bk.app, [*base, *extra, "--out", str(out)])
        assert res.exit_code == 2 and "refused" in res.output, (extra, res.output)
        assert nothing_written(out)
    ok = tmp_path / bk.R36_FILE
    res = runner.invoke(bk.app, [*base, "--form", "c2-s2", "--composition-from", "allowed-s1", "--captain-ruling",
                                 "R36", "--out", str(ok)])
    assert res.exit_code == 0, res.output
    assert "composition from allowed-s1 (row R36): the 5 scenario states at 0.3835" in res.output
    assert check_submission(ok, world.board).ok
    side = json.loads(ok.with_suffix(".json").read_text())
    assert side["combined_build"]["form"] == "c2-s2" and side["spec"]["composition"] == BLOCK
    assert side["combined_build"]["composition_line"]["captain_ruling"]["row"] == "R36"
    # the ruling row gone: refused
    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path, "", name="old.md"))
    res = runner.invoke(bk.app, [*base, "--form", "c2-s2", "--composition-from", "allowed-s1", "--captain-ruling",
                                 "R36", "--out", str(tmp_path / "n.h5ad")])
    assert res.exit_code == 2 and nothing_written(tmp_path / "n.h5ad")
