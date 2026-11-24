"""Row R22's third T3 search line in experiments/t3_decomp/build_two_knockout.py: the gata4_pooled weight g inside C2's
size-held product on the barred s2, `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`, only while
row R22's first cell names `--form c2-s2 --gata4-exponent 2.0`; one file, T3_c2_s2_g2.h5ad. The factor
phi_M ** c * phi_G ** (g c), c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), comes from a new function beside combine()
(g = 0 gives phi_M and g = 1 C2's factor, value for value); its guard (the file's den-weighted total within 5% of C2's
factor's, exit 6, nothing written) and its factor check (exit 5).

Offline, synthetic and fast: the throwaway world of tests/test_t3_two_knockout_build.py (a small synthetic WT file,
synthetic `transfer_mab21l2_all` and `gata4_pooled` tables, a decision file with the frozen hashes), a second run whose
gata4_pooled table is the world's at a quarter of its log size (it passes the guard at g = 2; the world's own fails
it), and small ledger files written under tmp_path. Nothing here reads the repository's own TRIES.md, builds from
released data or pins bytes built on one machine.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
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
PHI_M, PHI_G = th.PHI_M, th.PHI_G
QUARTER = PHI_G ** 0.25  # a milder gata4_pooled table: its file passes the guard at g = 2
LINE2 = ("by T3's second search line the T3 two-knockout builder is also to take `--gata4-exponent 2.0` on "
         "`--form c1-s2`, only with `--captain-ruling R22`")
LINE3 = ("by T3's third search line, on the captain's own yes naming T3_c2_s2_g2.h5ad, the T3 two-knockout builder is "
         "also to take `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`")
ROW = f"| R22 (the captain's order of a day: a search; ruling, no run, not a try; {LINE2}; {LINE3}) | T3 | x |\n"
ROW_LINE2 = f"| R22 (the captain's order of a day: a search; ruling, no run, not a try; {LINE2}) | T3 | x |\n"
OTHER_FORMS = ["identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c2-s1", "c2-s1x1.5",
               "c2-s1x2.0"]


def ledger(tmp_path: Path, rows: str = ROW, name: str = "TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text(f"# a throwaway ledger\n| R21 (another row) | T2 | y |\n{rows}| R41 (a ruling) | T3 | z |\n")
    return p


@pytest.fixture(scope="module")
def pw(world):
    """A run whose gata4_pooled table is QUARTER (same genes, same WT record): the line's file passes the guard."""
    run = world.d / "run_quarter"
    th.write_run(run, world.wt, {bk.TABLE_M: PHI_M, bk.TABLE_G: QUARTER})
    pins = {k: t3d.sha256(run / "factors" / f"{k}.csv") for k in (bk.TABLE_M, bk.TABLE_G)}
    return SimpleNamespace(run=run, pins=pins)


@pytest.fixture
def passing(pw, monkeypatch):
    """The quarter run's tables pinned (the world's are pinned by the autouse synthetic_pins)."""
    monkeypatch.setattr(bk, "TABLE_SHA256", dict(pw.pins))
    return pw


def build(world, out: Path, form: str = "c2-s2", **kw) -> dict:
    return bk.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=kw.pop("run_dir", world.run),
                         seed=SEED, threads=2, **kw)


def build_g2(world, tmp_path: Path, out: Path, **kw) -> dict:
    return build(world, out, gata4_exponent=2.0, captain_ruling="R22", tries=kw.pop("tries", ledger(tmp_path)), **kw)


def nothing_written(out: Path) -> bool:
    return not out.exists() and not out.with_suffix(".json").exists() and not Path(f"{out}.stamp.json").exists()


def loader_restored() -> bool:
    return bs.load_factor.__module__ == "build_submission" and bs.load_factor.__name__ == "load_factor"


class Spy:
    """Records the factor apply_gene_factor receives."""

    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


def den_of(world) -> np.ndarray:
    return np.asarray(bk._outside_form().wt_den(world.wt)[1], dtype=np.float64)


# ------------------------------------------------------------------------------------------------------ constants


def test_the_line_s_constants():
    assert (bk.C2G_FORM, bk.C2G_WEIGHTS, bk.C2G_OPTION, bk.C2G_FILE) == (
        "c2-s2", (2.0,), "`--form c2-s2 --gata4-exponent 2.0`", "T3_c2_s2_g2.h5ad")
    assert bk.C2G_MAX_DEN_SHIFT == 0.05 and bk.EXIT_GUARD_FAILED == 6 and bk.EXIT_FACTOR_MISMATCH == 5
    assert issubclass(bk.GuardFailed, bk.Refused) and not issubclass(bk.GuardFailed, bk.FactorMismatch)
    assert "R41" in bk.C2G_LIST_ENTRY and "R22" in bk.C2G_LIST_ENTRY
    # line 2's constants, the forms and combine's kinds are as they were: the line adds no form
    assert (bk.G2_ROW, bk.G2_EXPONENTS, bk.G2_FORM, bk.G2_OPTION, bk.G2_FILE) == (
        "R22", (2.0,), "c1-s2", "`--gata4-exponent 2.0`", "T3_c1_s2_g2.h5ad")
    assert bk.FORMS == ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2",
                        "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2")
    assert bk.KINDS == {"identity": 0.0, "c1": 1.0, "c2": 1.0}
    # line 2's option string is not inside line 3's, and line 3's is not line 2's: neither entry opens the other
    assert bk.G2_OPTION not in bk.C2G_OPTION and bk.C2G_OPTION not in LINE2
    assert not hasattr(bk, "captain_order_names") and not hasattr(bk, "wt_den")


# ------------------------------------------------------------------------------------------------------ the factor


def test_g_one_is_c2_s_factor_and_g_zero_phi_m_value_for_value():
    rng = np.random.default_rng(7)
    cases = [(PHI_M, PHI_G), (PHI_M, QUARTER), (np.array([2.0, 0.5, 1.0, 4.0]), np.array([2.0, 0.5, 1.0, 0.25]))]
    cases += [(np.exp(rng.normal(0, 0.6, 300)), np.exp(rng.normal(0.1, 0.4, 300))) for _ in range(5)]
    for m, g in cases:
        phi1, rec1 = bk.c2_gata4_weighted(m, g, 1.0)
        c2, rec_c2 = bk.combine("c2", 1.0, m, g)
        assert np.array_equal(phi1, c2) and phi1.tobytes() == c2.tobytes()  # value for value, bit for bit
        assert rec1["c"] == rec_c2["c"] and rec1["gata4_weight"] == rec_c2["gata4_weight"]
        assert rec1["mab21l2_exponent"] == rec_c2["mab21l2_exponent"]
        phi0, rec0 = bk.c2_gata4_weighted(m, g, 0.0)
        assert np.array_equal(phi0, np.asarray(m, dtype=np.float64)) and rec0["c"] == 1.0  # B's factor, phi_M
        assert rec0["gata4_weight"] == 0.0


def test_the_factor_by_hand_at_g_two():
    a = math.log(2.0)
    m = np.array([2.0, 0.5, 1.0, 4.0])
    g = np.array([2.0, 0.5, 1.0, 0.25])
    # ln m = a[1, -1, 0, 2], 2 ln g = a[2, -2, 0, -4]; sum = a[3, -3, 0, -2], mean -a/2
    sd_m = a * math.sqrt(5.0) / 2.0  # population SD of a[1, -1, 0, 2]
    sd_sum = a * math.sqrt((3.5 ** 2 + 2.5 ** 2 + 0.5 ** 2 + 1.5 ** 2) / 4.0)
    c = sd_m / sd_sum
    phi, rec = bk.c2_gata4_weighted(m, g, 2.0)
    assert rec["c"] == pytest.approx(c, rel=1e-14, abs=0) and rec["g"] == 2.0
    assert rec["mab21l2_exponent"] == pytest.approx(c, rel=1e-14) and rec["gata4_weight"] == pytest.approx(2 * c)
    assert np.allclose(phi, [2.0 ** (3 * c), 2.0 ** (-3 * c), 1.0, 2.0 ** (-2 * c)], rtol=1e-13, atol=0)
    assert np.allclose(phi, m ** c * g ** (2 * c), rtol=1e-15, atol=0)
    assert "weight g, not a power" in rec["c_rule"]


@pytest.mark.parametrize("g", [0.0, 0.5, 1.0, 2.0, 3.0])
def test_at_every_g_the_sd_of_ln_factor_is_phi_m_s_and_nothing_is_clipped(g):
    phi, rec = bk.c2_gata4_weighted(PHI_M, PHI_G, g)
    assert np.std(np.log(phi)) == pytest.approx(np.std(np.log(PHI_M)), rel=1e-12, abs=0)
    c = np.std(np.log(PHI_M)) / np.std(np.log(PHI_M) + g * np.log(PHI_G))
    assert rec["c"] == pytest.approx(c, rel=1e-13) and np.allclose(phi, np.exp(c * (np.log(PHI_M) + g * np.log(
        PHI_G))), rtol=1e-12, atol=0)
    if g == 0.0:
        assert phi.max() == 4.0 and phi.min() == 0.25  # phi_M's own clip ends, untouched


def test_the_factor_refuses_other_genes_bad_values_and_an_undefined_c():
    m = np.array([2.0, 0.5, 1.0])
    with pytest.raises(bk.Refused, match="same genes"):
        bk.c2_gata4_weighted(m, np.array([2.0, 0.5]), 2.0)
    with pytest.raises(bk.Refused, match="finite and positive"):
        bk.c2_gata4_weighted(m, np.array([2.0, 0.0, 1.0]), 2.0)
    with pytest.raises(bk.Refused, match="c is undefined"):
        bk.c2_gata4_weighted(np.ones(3), np.ones(3), 2.0)
    with pytest.raises(bk.Refused, match="not finite"):
        bk.c2_gata4_weighted(m, m, float("nan"))


def test_combine_is_unchanged():
    """combine on the world's tables, by hand (the existing tests' values, again): untouched by the line."""
    phi, rec = bk.combine("c1", 1.0, PHI_M, PHI_G)
    assert np.array_equal(phi, th.C1_HAND) and rec["gata4_weight"] == 1.0 and rec["c"] == 1.0
    lm, lg = np.log(PHI_M), np.log(PHI_G)
    for x in (1.0, 1.5, 2.0):
        phi, rec = bk.combine("c2", x, PHI_M, PHI_G)
        c = np.std(x * lm) / np.std(x * lm + lg)
        assert rec["c"] == pytest.approx(c, rel=1e-13) and np.array_equal(phi, (PHI_M ** (x * c)) * (PHI_G ** c))
    assert np.array_equal(bk.combine("identity", 1.0, PHI_M, PHI_G)[0], PHI_M)


# ------------------------------------------------------------------------------------------------------- the guard


def test_the_guard_by_hand():
    m = np.array([2.0, 0.5, 1.0, 1.0])
    g = np.array([2.0, 1.0, 4.0, 1.0])
    den = np.array([1.0, 2.0, 1.0, 4.0])
    r = bk.c2g_guard(m, g, 2.0, den)
    phi, rec = bk.c2_gata4_weighted(m, g, 2.0)
    c2, _ = bk.combine("c2", 1.0, m, g)
    s = den.sum()
    assert r["den_total_phi_m"] == pytest.approx(den @ m / s, rel=1e-14)
    assert r["den_total_c2"] == pytest.approx(den @ c2 / s, rel=1e-14)
    assert r["den_total_factor"] == pytest.approx(den @ phi / s, rel=1e-14)
    assert r["den_shift_vs_c2"] == pytest.approx((den @ phi) / (den @ c2) - 1.0, rel=1e-12)
    assert r["passed"] is (abs(r["den_shift_vs_c2"]) <= 0.05) and r["max_den_shift"] == 0.05
    assert r["den_shift_c2_vs_phi_m"] == pytest.approx((den @ c2) / (den @ m) - 1.0, rel=1e-12)
    assert r["den_shift_factor_vs_phi_m"] == pytest.approx((den @ phi) / (den @ m) - 1.0, rel=1e-12)
    assert r["vs_phi_m"] == "recorded; gates nothing" and r["gata4_weight_g"] == 2.0 and r["c"] == rec["c"]
    assert r["powers_applied"] == {bk.TABLE_M: rec["c"], bk.TABLE_G: 2.0 * rec["c"]}
    assert r["factor_sha256"] == hashlib.sha256(np.ascontiguousarray(phi, dtype=np.float64).tobytes()).hexdigest()
    with pytest.raises(bk.Refused, match="den"):
        bk.c2g_guard(m, g, 2.0, den[:-1])
    with pytest.raises(bk.Refused, match="den"):
        bk.c2g_guard(m, g, 2.0, -den)


def test_the_guard_passes_within_5_percent_of_c2_s_total_whatever_phi_m_s_shift_says(world):
    den = den_of(world)
    q = bk.c2g_guard(PHI_M, QUARTER, 2.0, den)  # the passing table: near C2, far from phi_M
    assert q["passed"] and abs(q["den_shift_vs_c2"]) <= 0.05 and q["den_shift_factor_vs_phi_m"] > 0.05
    w = bk.c2g_guard(PHI_M, PHI_G, 2.0, den)  # the world's own table: beyond 5% of C2's total
    assert not w["passed"] and abs(w["den_shift_vs_c2"]) > 0.05
    assert bk.c2g_guard(PHI_M, PHI_G, 1.0, den)["den_shift_vs_c2"] == 0.0  # C2 against itself
    # a constant gata4_pooled table: c = 1 at every g, the file is phi_M times 4 ** g, three times C2's total
    k = bk.c2g_guard(PHI_M, np.full(G, 4.0), 2.0, den)
    assert not k["passed"] and k["den_shift_vs_c2"] == pytest.approx(3.0, rel=1e-12)
    assert k["c"] == pytest.approx(1.0, rel=1e-14)


@pytest.mark.parametrize("shift,ok", [(0.0499, True), (-0.0499, True), (0.0, True), (0.0501, False),
                                      (-0.0501, False), (0.5, False)])
def test_the_guard_s_bound_is_5_percent_either_way(shift, ok, monkeypatch):
    """The bound on the shift against C2's total, both signs, at the edges (the totals stubbed)."""
    m, g, den = np.array([2.0, 0.5, 1.0]), np.array([1.5, 0.8, 1.1]), np.ones(3)
    real_combine = bk.combine

    def scaled(kind, x, pm, pg):
        phi, rec = real_combine(kind, x, pm, pg)
        file_phi, _ = bk.c2_gata4_weighted(pm, pg, 2.0)
        return file_phi * (np.sum(file_phi) / np.sum(file_phi * (1.0 + shift))), rec  # C2's total = file's / (1+s)

    monkeypatch.setattr(bk, "combine", scaled)
    r = bk.c2g_guard(m, g, 2.0, den)
    assert r["den_shift_vs_c2"] == pytest.approx(shift, abs=1e-12) and r["passed"] is ok


# ------------------------------------------------------------------------------------------------- the ruling row


@pytest.mark.parametrize("rows,ok", [
    (ROW, True),
    (ROW_LINE2, False),  # line 2's entry does not open line 3
    (ROW + ROW, False),  # two rows R22
    ("", False),  # no row R22
    (ROW.replace("| R22 (the captain's order", "| R22 (the captain's ruling"), False),
    (ROW.replace("ruling, no run", "a run"), False),
    (ROW.replace("`--form c2-s2 --gata4-exponent 2.0`", "--form c2-s2 --gata4-exponent 2.0"), False),  # no backticks
    (ROW.replace("`--form c2-s2 --gata4-exponent 2.0`", "`--form c2-s2 --gata4-exponent 3.0`"), False),
    (ROW.replace(f"; {LINE3}) | T3 | x |", ") | T3 | `--form c2-s2 --gata4-exponent 2.0` `--captain-ruling R22` |"),
     False),  # named outside the first cell only
    (ROW.replace("| R22 (", "| R22b ("), False),
])
def test_the_ruling_row_check(tmp_path, rows, ok):
    p = ledger(tmp_path, rows)
    if ok:
        rec = bk.check_c2_gata4("c2-s2", 2.0, "R22", p)
        assert rec == {"row": "R22", "option": "`--form c2-s2 --gata4-exponent 2.0`", "gata4_weight": 2.0,
                       "ledger": "experiments/harness/TRIES.md", "ledger_sha256": t3d.sha256(p)}
    else:
        with pytest.raises(bk.Refused, match="no single row R22"):
            bk.check_c2_gata4("c2-s2", 2.0, "R22", p)


def test_the_row_check_is_build_outside_form_s(tmp_path, monkeypatch):
    bo = bk._outside_form()
    seen = []
    monkeypatch.setattr(bo, "captain_order_names", lambda row, option, path=None: seen.append((row, option, path))
                        or True)
    p = tmp_path / "absent.md"
    monkeypatch.setattr(bs, "sha256", lambda path: "x")
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", p)["row"] == "R22"
    assert seen == [("R22", "`--form c2-s2 --gata4-exponent 2.0`", p)]


def test_line_3_s_entry_does_not_open_line_2(tmp_path):
    only3 = f"| R22 (the captain's order of a day: a search; ruling, no run, not a try; {LINE3}) | T3 | x |\n"
    with pytest.raises(bk.Refused, match="no single row R22"):
        bk.check_gata4("c1-s2", 2.0, "R22", ledger(tmp_path, only3))
    assert bk.check_gata4("c1-s2", 2.0, "R22", ledger(tmp_path, ROW))["option"] == "`--gata4-exponent 2.0`"


# ------------------------------------------------------------------------------------------------------ refusals


def test_without_the_option_on_c2_s2_nothing_of_the_line_is_read(monkeypatch):
    monkeypatch.setattr(bk, "_outside_form", lambda: pytest.fail("read on a path without the line"))
    for form in [*bk.FORMS, "c1-s1x1.5", "", None]:
        for kw in ({}, {"captain_ruling": "R22"}, {"captain_ruling": "R36"}):
            assert bk.check_c2_gata4(form, **kw) is None
    for form in [f for f in bk.FORMS if f != "c2-s2"]:
        assert bk.check_c2_gata4(form, 2.0, "R22") is None  # goes on to check_gata4, as before


@pytest.mark.parametrize("value", [3.0, 1.0, 0.0, 2.5, 1.5, 4.0, -2.0, float("nan"), float("inf"), "x", "", True,
                                   "2,0"])
def test_any_other_value_on_c2_s2_is_refused(tmp_path, value):
    with pytest.raises(bk.Refused, match="only 2.0"):
        bk.check_c2_gata4("c2-s2", value, "R22", ledger(tmp_path))


@pytest.mark.parametrize("ruling", ["R21", "R30", "R41", "r22", "R22 ", " R22", "R22b", "", "R36"])
def test_any_other_ruling_is_refused(tmp_path, ruling):
    with pytest.raises(bk.Refused, match="only row R22"):
        bk.check_c2_gata4("c2-s2", 2.0, ruling, ledger(tmp_path))


def test_the_option_without_the_ruling_is_refused(tmp_path):
    with pytest.raises(bk.Refused, match="only with --captain-ruling R22"):
        bk.check_c2_gata4("c2-s2", 2.0, None, ledger(tmp_path))


@pytest.mark.parametrize("form,kw,says", [
    ("c2-s2", {"gata4_exponent": 3.0, "captain_ruling": "R22"}, "only 2.0"),
    ("c2-s2", {"gata4_exponent": 2.0}, "only with --captain-ruling R22"),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R30"}, "only row R22"),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22", "rows": ROW_LINE2}, "no single row R22"),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22", "rows": ""}, "no single row R22"),
    ("c2-s2", {"captain_ruling": "R22"}, "taken only with --gata4-exponent 2.0 on --form c1-s2"),  # the ruling alone
    ("c2-s1", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "on --form c1-s2 only"),
    ("identity-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "on --form c1-s2 only"),
    ("c1-s1", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "on --form c1-s2 only"),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R36"}, "taken only with --composition-from"),
])
def test_the_builder_refuses_and_builds_nothing(world, tmp_path, monkeypatch, form, kw, says):
    for name in ("run_c2g_guard", "c2g_loader", "run_gata4_figures"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on a refusal"))
    kw = dict(kw)
    t = ledger(tmp_path, kw.pop("rows", ROW))
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match=says):
        build(world, out, form, tries=t, **kw)
    assert nothing_written(out) and loader_restored()


# --------------------------------------------------------------------------------------------- the guard at the build


def test_the_world_s_table_fails_the_guard_and_nothing_is_written(world, tmp_path, monkeypatch):
    monkeypatch.setattr(bk, "c2g_loader", lambda *a, **k: pytest.fail("the build was reached past a failed guard"))
    out = tmp_path / bk.C2G_FILE
    with pytest.raises(bk.GuardFailed, match="beyond 5%"):
        build_g2(world, tmp_path, out)
    assert nothing_written(out) and loader_restored()


def test_the_line_s_file_is_c2_s_command_with_the_weight_g_two(world, tmp_path, passing, monkeypatch):
    spy = Spy(monkeypatch)
    c2 = build(world, tmp_path / "c2.h5ad", run_dir=passing.run)
    out = tmp_path / bk.C2G_FILE
    side = build_g2(world, tmp_path, out, run_dir=passing.run)
    assert len(spy.seen) == 2
    want, rec = bk.c2_gata4_weighted(PHI_M, QUARTER, 2.0)
    assert np.array_equal(spy.seen[1], want) and np.array_equal(spy.seen[0], bk.combine("c2", 1.0, PHI_M, QUARTER)[0])
    # C2's command with one option pair: spec, inputs, streams, composition, seed, count and board are C2's
    assert side["spec"] == c2["spec"] and side["inputs"] == c2["inputs"] and side["spec_source"] == c2["spec_source"]
    assert side["inputs"]["rng_streams"] == ["submission:barred:s2", "submission:barred:s2:assemble"]
    assert side["construction"]["composition"] == c2["construction"]["composition"]
    assert side["construction"]["factor_exponent"] == 1.0 and side["output"]["check"]["ok"]
    cb, cb2 = side["combined_build"], c2["combined_build"]
    assert {k: v for k, v in cb.items() if k not in ("combination", "c2_gata4_line")} == {
        k: v for k, v in cb2.items() if k != "combination"}
    assert cb["form"] == "c2-s2" and cb["combination"]["g"] == 2.0 and cb["combination"]["c"] == rec["c"]
    line = cb["c2_gata4_line"]
    assert line["captain_ruling"]["row"] == "R22" and line["captain_ruling"]["option"] == bk.C2G_OPTION
    assert line["gata4_weight_g"] == 2.0 and line["file"] == "T3_c2_s2_g2.h5ad" and line["c"] == rec["c"]
    assert line["powers_applied"] == {bk.TABLE_M: rec["c"], bk.TABLE_G: 2.0 * rec["c"]}
    assert "weight g, not a power" in line["option_value"] and line["line"] == bk.C2G_LINE
    sha = hashlib.sha256(np.ascontiguousarray(want).tobytes()).hexdigest()
    assert line["factor_sha256_applied"] == line["guard"]["factor_sha256"] == sha
    # the guard's figures, by hand from den as build_outside_form.wt_den makes it
    g = line["guard"]
    den = den_of(world)
    s = den.sum()
    c2f = bk.combine("c2", 1.0, PHI_M, QUARTER)[0]
    assert np.array_equal(np.array(g["den"]), den) and "genes" not in g and g["n_genes"] == G
    assert g["den_total_c2"] == pytest.approx(den @ c2f / s, rel=1e-12)
    assert g["den_total_factor"] == pytest.approx(den @ want / s, rel=1e-12)
    assert g["den_shift_vs_c2"] == pytest.approx((den @ want) / (den @ c2f) - 1, rel=1e-12) and g["passed"]
    assert g["den_shift_factor_vs_phi_m"] > 0.05 and g["vs_phi_m"] == "recorded; gates nothing"  # it stops nothing
    assert g["tables_sha256"] == passing.pins and g["wt_sha256"] == t3d.sha256(world.wt)
    assert json.loads(out.with_suffix(".json").read_text())["combined_build"]["c2_gata4_line"]["file"] == bk.C2G_FILE
    assert json.loads(Path(f"{out}.stamp.json").read_text())["clean"] is True
    assert not any("cell" in k.lower() for k in line) and loader_restored()


def test_the_line_draws_c2_s_cells_and_differs_in_expression_and_builds_one_sha256(world, tmp_path, passing):
    import anndata as ad

    build(world, tmp_path / "c2.h5ad", run_dir=passing.run)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    sa = build_g2(world, tmp_path, tmp_path / "a" / bk.C2G_FILE, run_dir=passing.run)
    sb = build_g2(world, tmp_path, tmp_path / "b" / bk.C2G_FILE, run_dir=passing.run)
    assert (tmp_path / "a" / bk.C2G_FILE).read_bytes() == (tmp_path / "b" / bk.C2G_FILE).read_bytes()
    assert sa["output"]["sha256"] == sb["output"]["sha256"]
    a, c = ad.read_h5ad(tmp_path / "a" / bk.C2G_FILE), ad.read_h5ad(tmp_path / "c2.h5ad")
    assert list(a.obs_names) == list(c.obs_names) and a.obs.equals(c.obs)
    assert np.array_equal(a.obsm["spatial_3D"], c.obsm["spatial_3D"])
    assert not np.array_equal(np.asarray(a.X), np.asarray(c.X))


def test_a_factor_other_than_the_guard_s_exits_5_and_writes_nothing(world, tmp_path, passing, monkeypatch):
    real = bk.run_c2g_guard
    monkeypatch.setattr(bk, "run_c2g_guard", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.FactorMismatch, match="not the one the guard describes"):
        build_g2(world, tmp_path, out, run_dir=passing.run)
    assert nothing_written(out)
    monkeypatch.setattr(bk, "run_c2g_guard", lambda *a, **k: real(*a, **k) | {"genes": ["Z"]})
    with pytest.raises(bk.FactorMismatch, match="genes"):
        build_g2(world, tmp_path, out, run_dir=passing.run)
    assert nothing_written(out) and loader_restored()


def test_the_guard_reads_the_pinned_tables_only(world, tmp_path, passing, monkeypatch):
    monkeypatch.setattr(bk, "TABLE_SHA256", {**passing.pins, bk.TABLE_G: "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match="pinned table"):
        build_g2(world, tmp_path, out, run_dir=passing.run)
    assert nothing_written(out)


# ------------------------------------------------------------------------------- every existing path as it was


@pytest.mark.parametrize("form", list(bk.FORMS))
def test_every_existing_form_never_touches_the_line(world, tmp_path, monkeypatch, form):
    plain = build(world, tmp_path / "plain.h5ad", form)
    for name in ("check_c2_gata4", "c2_gata4_weighted", "c2g_guard", "run_c2g_guard", "c2g_loader"):
        if name != "check_c2_gata4":
            monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on {form}"))
    trapped = build(world, tmp_path / "trapped.h5ad", form, gata4_exponent=None, captain_ruling=None)
    assert (tmp_path / "plain.h5ad").read_bytes() == (tmp_path / "trapped.h5ad").read_bytes()
    assert plain["combined_build"] == trapped["combined_build"] and "c2_gata4_line" not in plain["combined_build"]


def test_line_2_s_path_is_unchanged_and_never_touches_line_3(world, tmp_path, monkeypatch):
    t = ledger(tmp_path)
    plain = build(world, tmp_path / "g2.h5ad", "c1-s2", gata4_exponent=2.0, captain_ruling="R22", tries=t)
    for name in ("c2_gata4_weighted", "c2g_guard", "run_c2g_guard", "c2g_loader"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on line 2"))
    trapped = build(world, tmp_path / "g2b.h5ad", "c1-s2", gata4_exponent=2.0, captain_ruling="R22", tries=t)
    assert (tmp_path / "g2.h5ad").read_bytes() == (tmp_path / "g2b.h5ad").read_bytes()
    gl = trapped["combined_build"]["gata4_line"]
    assert gl["factor_sha256_applied"] == plain["combined_build"]["gata4_line"]["factor_sha256_applied"]
    assert gl["captain_ruling"]["option"] == "`--gata4-exponent 2.0`" and "c2_gata4_line" not in trapped[
        "combined_build"]


# --------------------------------------------------------------------------------------------------------- CLI


def test_the_cli(world, tmp_path, passing, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    boards = th.cli_boards(tmp_path)

    def args(run, *extra, out):
        return ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
                str(boards), "--run-dir", str(run), "--seed", str(SEED), "--threads", "2", *extra, "--out", str(out)]

    runner = CliRunner()
    for extra in (["--form", "c2-s2", "--gata4-exponent", "3.0", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--gata4-exponent", "2.0"],
                  ["--form", "c2-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R41"],
                  ["--form", "c2-s1", "--gata4-exponent", "2.0", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--captain-ruling", "R22"]):
        out = tmp_path / "r.h5ad"
        res = runner.invoke(bk.app, args(passing.run, *extra, out=out))
        assert res.exit_code == 2 and "refused" in res.output, (extra, res.output)
        assert nothing_written(out)
    line = ["--form", "c2-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R22"]
    ok = tmp_path / bk.C2G_FILE
    res = runner.invoke(bk.app, args(passing.run, *line, out=ok))
    assert res.exit_code == 0, res.output
    assert "row R22's third T3 line: gata4_pooled weight g = 2 in C2's size-held product" in res.output
    assert json.loads(ok.with_suffix(".json").read_text())["combined_build"]["c2_gata4_line"]["gata4_weight_g"] == 2.0
    # the guard fails on the world's own tables: exit 6, nothing written
    monkeypatch.setattr(bk, "TABLE_SHA256", dict(world.pins))
    bad = tmp_path / "guard.h5ad"
    res = runner.invoke(bk.app, args(world.run, *line, out=bad))
    assert res.exit_code == 6 and "beyond 5%" in res.output and nothing_written(bad)
    # the factor is not the guard's: exit 5, nothing written
    monkeypatch.setattr(bk, "TABLE_SHA256", dict(passing.pins))
    real = bk.run_c2g_guard
    monkeypatch.setattr(bk, "run_c2g_guard", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    mis = tmp_path / "mismatch.h5ad"
    res = runner.invoke(bk.app, args(passing.run, *line, out=mis))
    assert res.exit_code == 5 and nothing_written(mis)
    monkeypatch.setattr(bk, "run_c2g_guard", real)
    # the ruling row naming line 2 only: refused
    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path, ROW_LINE2, name="line2.md"))
    res = runner.invoke(bk.app, args(passing.run, *line, out=tmp_path / "n.h5ad"))
    assert res.exit_code == 2 and nothing_written(tmp_path / "n.h5ad")
