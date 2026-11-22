"""Row R22's T3 search line in experiments/t3_decomp/build_outside_form.py: xall-s2 at --table-exponent 4.0, taken
only with --captain-ruling R22 while the ruling row's first cell names the option, after the line's guard (the
den-weighted total of phi_M * T ** 4 within 5% of phi_M's, from WT E8.75 and the pinned tables alone).

Synthetic data only: the WT file, phi_M, the tables, the listing and the ledger are throwaway fixtures (the helpers of
tests/test_t3_outside_forms.py and tests/test_t3_barred_exponent_build.py). Nothing here reads the repository's own
ledger for its state of today, and no built file's bytes are pinned: each comparison is between two builds made here."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
T3 = ROOT / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


to = sys.modules.get("_t3_outside_forms_helpers") or _load("_t3_outside_forms_helpers",
                                                           ROOT / "tests" / "test_t3_outside_forms.py")
bo, bs, ext, tx = to.bo, to.bs, to.ext, to.tx
t3d = bs.t3d
own_repo, lst, world, bw = to.own_repo, to.lst, to.world, to.bw  # fixtures, as the helpers define them
SEED, N, G, GENES, PHI = to.SEED, to.N, to.G, to.GENES, to.PHI

SENTENCE = ("by the lines added in decision log 2026-11-22 the T1 builder is also to take `--g-mom-damping 1.5` and "
            "`--g-mom-damping 2.0`, and the T3 outside-form builder `--table-exponent 4.0`, each only with "
            "`--captain-ruling R22`")
T1_PART = ("the builder is to take `--feng-fraction-mult 2.5` and `--feng-fraction-mult 3.0` only with "
           "`--captain-ruling R22`")
ROW = f"| R22 (the captain's order of a day: a search; ruling, no run, not a try; {T1_PART}; {SENTENCE}) | T1 | x |\n"


def ledger(tmp_path: Path, rows: str = ROW, name: str = "TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text(f"# a throwaway ledger\n| R21 (another row) | T2 | y |\n{rows}| K23 (a try) | T3 | z |\n")
    return p


def gentle() -> np.ndarray:
    """A table whose fourth power moves the den-weighted total by well under 5%: 12 genes at exp(+-0.01)."""
    T = np.ones(G)
    T[25:37] = np.exp(np.linspace(-0.01, 0.01, 12))
    return T


def strong() -> np.ndarray:
    """A table whose fourth power moves the den-weighted total far beyond 5%: 15 genes at 1.3."""
    T = np.ones(G)
    T[25:40] = 1.3
    return T


def tables_with(all_values) -> dict:
    return {**to.TABLES, "all": (all_values, ["GSE204908", "GSE277839"], ["unrelated"])}


def den_of(wt: Path) -> np.ndarray:
    from vec.prep import load_stage

    a = load_stage(wt, drop_duplicates=True)
    C = np.asarray(t3d.dense_counts(a), dtype=np.float64)
    return (C / np.maximum(C.sum(1, keepdims=True), 1)).mean(0)


def shift_of(wt: Path, T, e: float = 4.0) -> float:
    den = den_of(wt)
    return float((den @ (PHI * T ** e)) / (den @ PHI) - 1.0)


@pytest.fixture
def line(world, tmp_path, lst, monkeypatch):
    """phi_M's pin, the gentle table as `all` in a committed listing, a K23 reading naming x_e2 and a ledger whose
    row R22 names the option."""
    monkeypatch.setattr(bo, "TABLE_SHA256", world.pin)
    s = to.setup_tables(tmp_path, lst, monkeypatch, tables=tables_with(gentle()))
    s.k23 = to.k23_reading(lst, names=("x_e2",))
    s.tries = ledger(tmp_path)
    return s


def build(world, s, out: Path, **kw) -> dict:
    return bo.build_form(world.dec, kw.pop("form", "xall-s2"), world.wt, world.board, N, out, run_dir=world.run,
                         seed=SEED, threads=2, env=s.env, k23_reading=kw.pop("k23_reading", s.k23),
                         tries=kw.pop("tries", s.tries), **kw)


def frozen_with(world, tmp_path, monkeypatch, factor, slot="s2", name="frozen.h5ad") -> bytes:
    original = bs.load_factor
    with monkeypatch.context() as m:
        m.setattr(bs, "load_factor", lambda run_dir, nm, genes, dec: (lambda q: (q[0] * factor, q[1]))(
            original(run_dir, nm, genes, dec)))
        bs.build_submission(world.dec, "barred", slot, world.wt, world.board, N, tmp_path / name, run_dir=world.run,
                            seed=SEED, threads=2)
    return (tmp_path / name).read_bytes()


# ------------------------------------------------------------------------------------------------ constants


def test_the_line_s_constants():
    assert (bo.R22_ROW, bo.R22_EXPONENT, bo.R22_FORM, bo.R22_OPTION) == ("R22", 4.0, "xall-s2",
                                                                         "`--table-exponent 4.0`")
    assert bo.GUARD_MAX_DEN_SHIFT == 0.05 and bo.GUARD_TABLE == "all" and bo.EXIT_GUARD_REFUSED == 5
    assert bo.BASES["s2"] == ("s2", 1.0, "1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f")
    assert bo.TABLE_SHA256.startswith("64b1e0b8572c")  # patched only inside fixtures
    assert bo.TRIES_FILE == ROOT / "experiments" / "harness" / "TRIES.md"
    assert issubclass(bo.GuardRefused, bo.Refused) and not issubclass(bo.GuardRefused, bo.FormRefused)
    assert bo.TABLE_EXPONENTS == (1.0, 2.0, 4.0)


def test_the_guard_by_hand():
    den = np.linspace(1.0, 2.0, G)
    T = strong()
    g = bo.exponent_guard(PHI, 1.0, T, 4.0, den)
    want = float((den @ (PHI * T ** 4)) / (den @ PHI) - 1.0)
    assert np.isclose(g["den_shift"], want, rtol=1e-12) and not g["passed"] and g["max_den_shift"] == 0.05
    g1 = bo.exponent_guard(PHI, 1.0, gentle(), 4.0, den)
    assert g1["passed"] and abs(g1["den_shift"]) < 0.05
    # exactly at the bound passes, just beyond fails (the rule: within 5%)
    one = np.ones(G)
    edge = np.ones(G)
    edge[:] = 1.0499 ** 0.25
    assert bo.exponent_guard(one, 1.0, edge, 4.0, np.ones(G))["passed"]
    edge[:] = 1.0501 ** 0.25
    assert not bo.exponent_guard(one, 1.0, edge, 4.0, np.ones(G))["passed"]
    edge[:] = 0.9499 ** 0.25
    assert not bo.exponent_guard(one, 1.0, edge, 4.0, np.ones(G))["passed"]


# ------------------------------------------------------------------------------------------------ the ruling row


@pytest.mark.parametrize("rows, ok", [
    (ROW, True),
    (ROW.replace(f"; {SENTENCE}", ""), False),  # the row before the sentence is added: T1's options only
    (ROW + ROW, False),  # two rows R22
    ("", False),  # no row R22
    (ROW.replace("| R22 (the captain's order", "| R22 (the captain's ruling"), False),
    (ROW.replace("ruling, no run", "a run"), False),
    (ROW.replace(SENTENCE, SENTENCE.replace("`--table-exponent 4.0`", "--table-exponent 4.0")).replace(
        T1_PART, "the builder"), False),  # the option without backticks
    (ROW.replace(f"; {T1_PART}; {SENTENCE}) | T1 | x |",
                 ") | T1 | `--table-exponent 4.0` with `--captain-ruling R22` |"),
     False),  # the option only in a later cell
    (ROW.replace(T1_PART, "the T1 values").replace("each only with `--captain-ruling R22`",
                                                   "each only with a flag"),
     False),  # no `--captain-ruling R22` in the first cell
    (ROW.replace("| R22 (", "| R22b ("), False),
])
def test_the_ruling_row_check(tmp_path, rows, ok):
    p = ledger(tmp_path, rows)
    assert bo.captain_order_names("R22", "`--table-exponent 4.0`", p) is ok
    if ok:
        assert bo.check_captain_ruling("R22", "xall-s2", 4.0, p)["ledger_sha256"] == t3d.sha256(p)
    else:
        with pytest.raises(bo.Refused, match="no single row R22"):
            bo.check_captain_ruling("R22", "xall-s2", 4.0, p)
    assert bo.captain_order_names("R22", "`--table-exponent 4.0`", tmp_path / "absent.md") is False


def test_the_flag_takes_r22_xall_s2_at_4_only(tmp_path):
    p = ledger(tmp_path)
    for ruling in ("R20", "R21", "R23", "r22", "", "K23"):
        with pytest.raises(bo.Refused, match="only row R22"):
            bo.check_captain_ruling(ruling, "xall-s2", 4.0, p)
    for form, e in (("xall-s2", 2.0), ("xall-s2", 1.0), ("xall-s1", 4.0), ("xall-s1x2.0", 4.0), ("xunrel-s2", 4.0),
                    ("identity-s2", 4.0)):
        with pytest.raises(bo.Refused, match="takes --form xall-s2 --table-exponent 4 only"):
            bo.check_captain_ruling("R22", form, e, p)


# ------------------------------------------------------------------------------------------------ the builder


def test_4_is_refused_without_the_flag(world, line, tmp_path):
    """Without --captain-ruling the K23 path holds as before: 4.0 needs a K23 reading naming x_e4 (this one names
    x_e2), whatever the ledger says."""
    with pytest.raises(bo.Refused, match="does not name x_e4"):
        build(world, line, tmp_path / "x.h5ad", table_exponent=4.0)
    assert not (tmp_path / "x.h5ad").exists() and not (tmp_path / "x.json").exists()


def test_4_is_accepted_with_the_flag(world, line, tmp_path, monkeypatch):
    spy = to.Spy(monkeypatch)
    side = build(world, line, tmp_path / "e4.h5ad", table_exponent=4.0, captain_ruling="R22")
    (phi,) = spy.seen
    T = gentle()
    assert np.array_equal(phi, ext.product(np.asarray(PHI, dtype=np.float64), 1.0, T, 4.0))
    ob = side["outside_build"]
    assert ob["form"] == "xall-s2" and ob["exponent"] == 4.0 and ob["base"]["x"] == 1.0
    assert ob["base"]["reference_sha256"] == bo.BASES["s2"][2]
    assert ob["captain_ruling"]["row"] == "R22" and ob["captain_ruling"]["option"] == "`--table-exponent 4.0`"
    assert ob["captain_ruling"]["ledger_sha256"] == t3d.sha256(line.tries)
    assert ob["line"] == bo.R22_LINE and ob["list_entry"] == bo.R22_LIST_ENTRY
    assert ob["priced_at"] == "x = 1 only" and "exponent 1 only" in ob["priced_note"]
    assert "row K23 did not name exponent 4" in ob["priced_note"]
    assert ob["k23_reading"] is None and ob["k22_reading"] is None and ob["table"]["sha256"] == line.shas["all"]
    g = ob["guard"]
    assert g["passed"] and g["max_den_shift"] == 0.05 and "genes" not in g
    assert np.isclose(g["den_shift"], shift_of(world.wt, T), rtol=1e-9, atol=1e-15)
    assert g["wt_sha256"] == t3d.sha256(world.wt) and g["table_m_sha256"] == world.pin
    assert g["table_sha256"] == line.shas["all"] and g["n_genes"] == G
    assert json.loads((tmp_path / "e4.json").read_text())["outside_build"]["guard"] == g
    assert side["construction"]["factor_exponent"] == 1.0
    # the file is the frozen builder's barred s2 given phi_M * T ** 4, and a second build gives the same bytes
    assert frozen_with(world, tmp_path, monkeypatch, T ** 4.0) == (tmp_path / "e4.h5ad").read_bytes()
    build(world, line, tmp_path / "e4b.h5ad", table_exponent=4.0, captain_ruling="R22")
    assert (tmp_path / "e4b.h5ad").read_bytes() == (tmp_path / "e4.h5ad").read_bytes()


def test_the_flag_is_refused_without_the_row(world, line, tmp_path):
    for rows in (ROW.replace(f"; {SENTENCE}", ""), ""):
        t = ledger(tmp_path, rows, name="other.md")
        with pytest.raises(bo.Refused, match="no single row R22"):
            build(world, line, tmp_path / "x.h5ad", table_exponent=4.0, captain_ruling="R22", tries=t)
    assert not (tmp_path / "x.h5ad").exists() and not (tmp_path / "x.json").exists()


def test_the_guard_failing_refuses_and_writes_nothing(world, tmp_path, lst, monkeypatch):
    monkeypatch.setattr(bo, "TABLE_SHA256", world.pin)
    s = to.setup_tables(tmp_path, lst, monkeypatch, tables=tables_with(strong()))
    s.k23, s.tries = to.k23_reading(lst, names=("x_e2",)), ledger(tmp_path)
    assert shift_of(world.wt, strong()) > 0.05  # the fixture really is beyond the bound

    def never(*a, **k):
        raise AssertionError("the frozen builder ran after a failed guard")

    monkeypatch.setattr(bs, "build_submission", never)
    with pytest.raises(bo.GuardRefused, match="beyond 5%"):
        build(world, s, tmp_path / "x.h5ad", table_exponent=4.0, captain_ruling="R22")
    assert not list(tmp_path.glob("x.*"))


def test_other_values_stay_refused(world, line, tmp_path):
    for e in (0.5, 1.5, 3.0, 8.0, -4.0, "x"):
        for ruling in (None, "R22"):
            with pytest.raises(bo.Refused, match="--table-exponent"):
                build(world, line, tmp_path / "x.h5ad", table_exponent=e, captain_ruling=ruling)
    for e in (1.0, 2.0):
        with pytest.raises(bo.Refused, match="takes --form xall-s2 --table-exponent 4 only"):
            build(world, line, tmp_path / "x.h5ad", table_exponent=e, captain_ruling="R22")
    for form in ("xall-s1", "xall-s1x1.5", "xall-s1x2.0"):
        with pytest.raises(bo.Refused, match="takes --form xall-s2"):
            build(world, line, tmp_path / "x.h5ad", form=form, table_exponent=4.0, captain_ruling="R22")
    for form in ("xunrel-s2", "identity-s2"):
        with pytest.raises(bo.Refused, match="xall only"):
            build(world, line, tmp_path / "x.h5ad", form=form, table_exponent=4.0, captain_ruling="R22")
    with pytest.raises(bo.Refused, match="only row R22"):
        build(world, line, tmp_path / "x.h5ad", table_exponent=4.0, captain_ruling="R20")
    assert not list(tmp_path.glob("x.*"))


def test_the_existing_exponent_s_file_reproduces(world, line, tmp_path, monkeypatch):
    """The queued file's command (xall-s2 at 2.0, no flag, the K23 reading naming x_e2) and the base's identity build
    as before on the fixture: the frozen builder's barred s2 given phi_M * T ** 2 and given nothing, byte for byte,
    and the sidecar block is row K23's (no key of the R22 line)."""
    T = gentle()
    side = build(world, line, tmp_path / "e2.h5ad", table_exponent=2.0)
    assert frozen_with(world, tmp_path, monkeypatch, T ** 2.0) == (tmp_path / "e2.h5ad").read_bytes()
    ob = side["outside_build"]
    assert ob["exponent"] == 2.0 and ob["k23_reading"]["names"] == ["x_e2"] and ob["list_entry"] == bo.K23_LIST_ENTRY
    assert set(ob) == to.OUTSIDE_BUILD_KEYS | {"k23_reading", "priced_at"}
    build(world, line, tmp_path / "id.h5ad", form="identity-s2")
    frozen = bs.build_submission(world.dec, "barred", "s2", world.wt, world.board, N, tmp_path / "s2.h5ad",
                                 run_dir=world.run, seed=SEED, threads=2)
    assert frozen and (tmp_path / "s2.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()


# ------------------------------------------------------------------------------------------------ the command line


def test_the_cli_exit_codes(world, line, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bo, "READING_K23", line.k23)
    monkeypatch.setattr(bo, "TRIES_FILE", line.tries)
    monkeypatch.setenv(ext.ENV, line.env[ext.ENV])
    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(to._boards(tmp_path)), "--run-dir", str(world.run), "--threads", "2", "--form", "xall-s2",
            "--table-exponent", "4.0"]
    runner = CliRunner()
    res = runner.invoke(bo.app, [*base, "--out", str(tmp_path / "no.h5ad")])
    assert res.exit_code == 2 and "does not name x_e4" in res.output and not (tmp_path / "no.h5ad").exists()
    res = runner.invoke(bo.app, [*base, "--captain-ruling", "R20", "--out", str(tmp_path / "r20.h5ad")])
    assert res.exit_code == 2 and not (tmp_path / "r20.h5ad").exists()
    res = runner.invoke(bo.app, [*base, "--captain-ruling", "R22", "--out", str(tmp_path / "ok.h5ad")])
    assert res.exit_code == 0, res.output
    assert json.loads((tmp_path / "ok.json").read_text())["outside_build"]["captain_ruling"]["row"] == "R22"
    monkeypatch.setattr(bo, "TRIES_FILE", ledger(tmp_path, ROW.replace(f"; {SENTENCE}", ""), name="old.md"))
    res = runner.invoke(bo.app, [*base, "--captain-ruling", "R22", "--out", str(tmp_path / "old.h5ad")])
    assert res.exit_code == 2 and not (tmp_path / "old.h5ad").exists()


def test_the_cli_s_guard_exit_code(world, tmp_path, lst, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bo, "TABLE_SHA256", world.pin)
    s = to.setup_tables(tmp_path, lst, monkeypatch, tables=tables_with(strong()))
    monkeypatch.setattr(bo, "TRIES_FILE", ledger(tmp_path))
    monkeypatch.setenv(ext.ENV, s.env[ext.ENV])
    res = CliRunner().invoke(bo.app, ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic",
                                      "--boards-dir", str(to._boards(tmp_path)), "--run-dir", str(world.run),
                                      "--threads", "2", "--form", "xall-s2", "--table-exponent", "4.0",
                                      "--captain-ruling", "R22", "--out", str(tmp_path / "g.h5ad")])
    assert res.exit_code == 5 and "the line's guard" in res.output
    assert not (tmp_path / "g.h5ad").exists() and not (tmp_path / "g.json").exists()
