"""The options on C2's command in experiments/t3_decomp/build_two_knockout.py (the T3 design of the design round of
2026-11-24, section 2, as its cold challenge fixes it): on `--form c2-s2`, `--gata4-exponent inf` (F1,
T3_c2_s2_gonly.h5ad: the gata4_pooled table alone at C2's held size, phi_M dropped), `--size-scale 0.67` (F2,
T3_c2_s2_s067.h5ad: C2's ln factor times 0.67) and both (F3, T3_c2_s2_gonly_s067.h5ad), each only with
`--captain-ruling R22` and only while row R22's first cell holds the file's option string as one whole
backtick-delimited string: `--form c2-s2 --gata4-exponent inf`, `--form c2-s2 --size-scale 0.67`,
`--form c2-s2 --gata4-exponent inf --size-scale 0.67`. The factor against its formula; the range guard and line 3's 5%
den check (exit 6, nothing written), the latter lifted only by a row named with `--den-ruling` that names the file's
string and says "the 5% check gates nothing for it"; exit 5; every refusal; and every other command untouched. No
cosine or correlation of the two tables is computed here.

Offline, synthetic and fast: the throwaway world of tests/test_t3_two_knockout_build.py (a small synthetic WT file,
synthetic `transfer_mab21l2_all` and `gata4_pooled` tables, a decision file with the frozen hashes), runs whose two
tables are built here so that each check passes or fails by construction, and small ledger files written under
tmp_path. Nothing here reads the repository's own TRIES.md, builds from released data or pins bytes built on one
machine.
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
INF = math.inf


def alternating(hi: float, lo: float, start: int, stop: int) -> np.ndarray:
    """1 on every gene but start..stop-1, where it alternates hi, lo, hi, lo."""
    v = np.ones(G)
    v[start:stop:2], v[start + 1:stop:2] = hi, lo
    return v


# Tables whose verdicts follow by hand. MILD: every file's factor and C2's lie within [0.98, 1.02], so both checks pass
# for every file at any den. ORTH: phi_M moves genes 0-9 and phi_G genes 10-19 by x2 / x0.5, so c_inf = 1 (F1 is phi_G)
# and C2 = exp((ln phi_M + ln phi_G) / sqrt 2): F1 spans [0.5, 2] against C2's [2 ** -0.707, 2 ** 0.707]. WIDE: the same
# at x16, where F1's [1/16, 16] passes twice C2's 16 ** 0.707. SPIKE: phi_G moves one gene, c_inf ~ 3.2, so F1 and F3
# are too wide. SHIFT: phi_G is 4 on all genes but one, c_inf ~ 1.6: each file's factor sits far from C2's gene by gene
# (F1 and F3 above it, F2 below it), so each file's den-weighted total is far from C2's at any den; the range passes.
MILD = (alternating(1.02, 1 / 1.02, 0, 10), alternating(1.02, 1 / 1.02, 10, 20))
ORTH = (alternating(2.0, 0.5, 0, 10), alternating(2.0, 0.5, 10, 20))
WIDE = (alternating(16.0, 0.0625, 0, 10), alternating(16.0, 0.0625, 10, 20))
SPIKE = (alternating(4.0, 0.25, 0, 10), np.where(np.arange(G) == 10, 4.0, 1.0))
SHIFT = (alternating(2.0, 0.5, 0, 10), np.where(np.arange(G) == 39, 1.0, 4.0))
TABLES = {"mild": MILD, "orth": ORTH, "wide": WIDE, "spike": SPIKE, "shift": SHIFT}
RANGE_OK = {"mild": {"F1": True, "F2": True, "F3": True}, "orth": {"F1": True, "F2": True, "F3": True},
            "wide": {"F1": False, "F2": True, "F3": True}, "spike": {"F1": False, "F2": True, "F3": False},
            "shift": {"F1": True, "F2": True, "F3": True}}

OPT = {"F1": "--form c2-s2 --gata4-exponent inf", "F2": "--form c2-s2 --size-scale 0.67",
       "F3": "--form c2-s2 --gata4-exponent inf --size-scale 0.67"}
LINE2 = ("by T3's second search line the T3 two-knockout builder is also to take `--gata4-exponent 2.0` on "
         "`--form c1-s2`, only with `--captain-ruling R22`")
LINE3 = ("by T3's third search line the T3 two-knockout builder is also to take `--form c2-s2 --gata4-exponent 2.0`, "
         "only with `--captain-ruling R22`")


def new_entry(*names: str) -> str:
    return ("the T3 two-knockout builder is also to take " + ", ".join(f"`{OPT[n]}` ({n})" for n in names)
            + ", each only with `--captain-ruling R22`")


OPENER = "| R22 (the captain's order of a day: a search; ruling, no run, not a try"
ROW = f"{OPENER}; {LINE2}; {LINE3}; {new_entry('F1', 'F2', 'F3')}) | T3 | x |\n"
ROW_OLD = f"{OPENER}; {LINE2}; {LINE3}) | T3 | x |\n"  # the old lines only
ROW_NEW = f"{OPENER}; {new_entry('F1', 'F2', 'F3')}) | T3 | x |\n"  # the new strings only


def row_only(*names: str) -> str:
    return f"{OPENER}; {new_entry(*names)}) | T3 | x |\n"


def lift_row(*names: str, row: str = "R44", words: str = "the 5% check gates nothing for it",
             ruling: str = "ruling, no run, not a try") -> str:
    files = " and ".join(f"`{OPT[n]}`" for n in names)
    return f"| {row} (the captain's yes of a day; {ruling}; for {files}: {words}) | T3 | y |\n"


OPTS = {"F1": {"gata4_exponent": INF}, "F2": {"size_scale": 0.67}, "F3": {"gata4_exponent": INF, "size_scale": 0.67}}
FILES = {"F1": "T3_c2_s2_gonly.h5ad", "F2": "T3_c2_s2_s067.h5ad", "F3": "T3_c2_s2_gonly_s067.h5ad"}
POINTS = {"F1": (INF, 1.0), "F2": (1.0, 0.67), "F3": (INF, 0.67)}
LINE_FUNCS = ("c2_gonly_factor", "gonly_guards", "run_gonly_guards", "gonly_loader", "gonly_failure", "check_den_lift")


def ledger(tmp_path: Path, rows: str = ROW, name: str = "TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text(f"# a throwaway ledger\n| R21 (another row) | T2 | y |\n{rows}| R41 (a ruling) | T3 | z |\n")
    return p


@pytest.fixture(scope="module")
def runs(world):
    """One run per table set (same genes, same WT record) and its two tables' sha256."""
    out = {}
    for key, (m, g) in TABLES.items():
        run = world.d / f"run_{key}"
        th.write_run(run, world.wt, {bk.TABLE_M: m, bk.TABLE_G: g})
        pins = {k: t3d.sha256(run / "factors" / f"{k}.csv") for k in (bk.TABLE_M, bk.TABLE_G)}
        out[key] = SimpleNamespace(run=run, pins=pins, m=m, g=g)
    return out


def use(runs, monkeypatch, key: str):
    """The run `key`'s tables pinned (the world's are pinned by the autouse synthetic_pins)."""
    monkeypatch.setattr(bk, "TABLE_SHA256", dict(runs[key].pins))
    return runs[key]


def build(world, out: Path, form: str = "c2-s2", **kw) -> dict:
    return bk.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=kw.pop("run_dir", world.run),
                         seed=SEED, threads=2, **kw)


def build_f(world, tmp_path: Path, name: str, out: Path, rows: str = ROW, **kw) -> dict:
    return build(world, out, captain_ruling="R22", tries=ledger(tmp_path, rows, name=f"{out.stem}_TRIES.md"),
                 **OPTS[name], **kw)


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


def sd(v) -> float:
    v = np.asarray(v, dtype=np.float64)
    return float(np.sqrt(np.mean((v - v.mean()) ** 2)))


def ln_factor_by_hand(m, g_tab, g: float, s: float = 1.0) -> np.ndarray:
    """The family's ln factor, written out: s c_g (ln phi_M + g ln phi_G); at g = inf s (SD(ln M) / SD(ln G)) ln G."""
    lm, lg = np.log(m), np.log(g_tab)
    if g == INF:
        return s * (sd(lm) / sd(lg)) * lg
    return s * (sd(lm) / sd(lm + g * lg)) * (lm + g * lg)


# ------------------------------------------------------------------------------------------------------ constants


def test_the_line_s_constants():
    assert (bk.GONLY_FORM, bk.GONLY_SCALES, bk.GONLY_OPTIONS, bk.GONLY_FILES, bk.GONLY_POINTS) == (
        "c2-s2", (0.67,), OPT, FILES, POINTS)
    assert bk.GONLY_RANGE_FACTOR == 2.0 and bk.GONLY_LIFT_WORDS == "the 5% check gates nothing for it"
    assert bk.C2G_MAX_DEN_SHIFT == 0.05  # line 3's 5% check, the one these files get
    assert bk.EXIT_GUARD_FAILED == 6 and bk.EXIT_FACTOR_MISMATCH == 5 and bk.EXIT_REFUSED == 2
    assert issubclass(bk.GonlyGuardFailed, bk.GuardFailed) and not issubclass(bk.GonlyGuardFailed, bk.FactorMismatch)
    # no distinctness guard: nothing in the builder computes a cosine of one factor with another
    for gone in ("GONLY_MAX_COSINE", "GONLY_REFERENCES", "GONLY_DISTINCT_RULE", "GONLY_NOT_CHECKED", "_cosine",
                 "GONLY_OPTION_G", "GONLY_OPTION_S", "GONLY_NAMES"):
        assert not hasattr(bk, gone), gone
    src = (ROOT / "experiments" / "t3_decomp" / "build_two_knockout.py").read_text()
    line_src = src[src.index("def _is_pos_inf"):src.index("def build_form")]
    assert "cos" not in line_src.lower().replace("cost", "") and "corrcoef" not in line_src
    # the earlier lines' constants, the forms and combine's kinds are as they were: the line adds no form
    assert (bk.G2_OPTION, bk.C2G_OPTION, bk.C2G_WEIGHTS, bk.G2_EXPONENTS) == (
        "`--gata4-exponent 2.0`", "`--form c2-s2 --gata4-exponent 2.0`", (2.0,), (2.0,))
    assert bk.FORMS == ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2",
                        "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2")
    assert bk.KINDS == {"identity": 0.0, "c1": 1.0, "c2": 1.0}
    assert not hasattr(bk, "captain_order_names") and not hasattr(bk, "wt_den")


def test_the_whole_string_match_and_the_containment_hazard():
    """F1's string sits inside F3's: containment alone could open F1 on F3's entry; the whole-string match does not.
    No new string is or contains line 2's or line 3's, and none of those is inside a new one."""
    assert OPT["F1"] in OPT["F3"]
    cell_f3 = f"R22 (...; ruling, no run; `{OPT['F3']}` and `--captain-ruling R22`)"
    assert bk.names_whole(cell_f3, OPT["F3"]) and not bk.names_whole(cell_f3, OPT["F1"])
    cell_f1 = f"R22 (...; `{OPT['F1']}`)"
    assert bk.names_whole(cell_f1, OPT["F1"]) and not bk.names_whole(cell_f1, OPT["F3"])
    odd = f"R22 (`x`{OPT['F1']}`)"  # the backticked string is in the text, but not as one whole delimited string
    assert f"`{OPT['F1']}`" in odd and not bk.names_whole(odd, OPT["F1"])
    assert not bk.names_whole(None, OPT["F1"]) and not bk.names_whole(f"`{OPT['F1']} `", OPT["F1"])
    old = [bk.G2_OPTION.strip("`"), bk.C2G_OPTION.strip("`"), bk._outside_form().R22_OPTION.strip("`"),
           *(n.strip("`") for n in bk.R36_NAMES)]
    for a in OPT.values():
        for b in old:
            assert a != b and f"`{a}`" not in f"`{b}`" and f"`{b}`" not in f"`{a}`", (a, b)


# ------------------------------------------------------------------------------------------------------ the factor


@pytest.mark.parametrize("seed", range(6))
def test_the_factor_at_g_inf_is_the_formula(seed):
    rng = np.random.default_rng(seed)
    m = np.exp(rng.normal(0.05, 0.6, 300))
    g = np.exp(rng.normal(-0.02, 0.35, 300))
    phi, rec = bk.c2_gonly_factor(m, g, INF)
    c_inf = sd(np.log(m)) / sd(np.log(g))
    assert rec["c"] == pytest.approx(c_inf, rel=1e-12) and rec["g"] == "inf" and rec["size_scale"] == 1.0
    assert np.allclose(np.log(phi), c_inf * np.log(g), rtol=1e-12, atol=1e-15)
    assert np.allclose(np.log(phi), ln_factor_by_hand(m, g, INF), rtol=1e-12, atol=1e-15)
    assert np.array_equal(phi, g ** rec["c"])  # the operation itself, bit for bit
    assert sd(np.log(phi)) == pytest.approx(sd(np.log(m)), rel=1e-12)  # held at C2's size: phi_M's SD
    assert rec["mab21l2_exponent"] == 0.0 and rec["gata4_weight"] == pytest.approx(c_inf, rel=1e-12)
    # phi_M enters by its SD only (dropped otherwise): a permutation of it gives the same factor
    assert np.array_equal(bk.c2_gonly_factor(rng.permutation(m), g, INF)[0], phi)
    # the family's limit: the size-held product at a large finite g approaches it
    far, _ = bk.c2_gata4_weighted(m, g, 1e7)
    assert np.allclose(np.log(far), np.log(phi), rtol=1e-5, atol=1e-6)
    assert "corr_ln_m_ln_g" not in rec


def test_the_factor_at_g_inf_by_hand():
    m = np.array([2.0, 0.5, 1.0, 4.0])  # ln = a[1, -1, 0, 2]: SD a sqrt(5) / 2
    g = np.array([2.0, 0.5, 2.0, 0.5])  # ln = a[1, -1, 1, -1]: SD a
    c_inf = math.sqrt(5.0) / 2.0
    phi, rec = bk.c2_gonly_factor(m, g, INF)
    assert rec["c"] == pytest.approx(c_inf, rel=1e-14)
    assert np.allclose(phi, [2.0 ** c_inf, 2.0 ** -c_inf, 2.0 ** c_inf, 2.0 ** -c_inf], rtol=1e-14, atol=0)
    assert "phi_M dropped" in rec["c_rule"]


@pytest.mark.parametrize("g", [INF, 1.0])
def test_the_size_scale_multiplies_the_ln_factor_after_the_size_held_combination(g):
    rng = np.random.default_rng(11)
    for m, t in [(PHI_M, PHI_G), ORTH, SHIFT, (np.exp(rng.normal(0, 0.6, 200)), np.exp(rng.normal(0, 0.4, 200)))]:
        held, rec1 = bk.c2_gonly_factor(m, t, g)
        phi, rec = bk.c2_gonly_factor(m, t, g, 0.67)
        assert np.allclose(np.log(phi), 0.67 * np.log(held), rtol=1e-12, atol=1e-15)
        assert np.allclose(np.log(phi), ln_factor_by_hand(m, t, g, 0.67), rtol=1e-11, atol=1e-14)
        assert np.array_equal(phi, held ** 0.67)  # after the combination, bit for bit
        assert rec["size_scale"] == 0.67 and rec["c"] == rec1["c"]
        assert rec["gata4_weight"] == pytest.approx(0.67 * rec1["gata4_weight"], rel=1e-15)
        assert rec["mab21l2_exponent"] == pytest.approx(0.67 * rec1["mab21l2_exponent"], rel=1e-15)
        assert sd(np.log(phi)) == pytest.approx(0.67 * sd(np.log(m)), rel=1e-12)


def test_g_one_without_a_size_scale_is_c2_s_factor_value_for_value():
    for m, t in [(PHI_M, PHI_G), ORTH, WIDE, SPIKE, SHIFT]:
        phi, rec = bk.c2_gonly_factor(m, t, 1.0)
        c2, rc2 = bk.combine("c2", 1.0, m, t)
        assert phi.tobytes() == c2.tobytes() and rec["c"] == rc2["c"]
        assert rec["mab21l2_exponent"] == rec["gata4_weight"] == rc2["c"]
        assert bk.c2_gonly_factor(m, t, 1.0, 0.67)[0].tobytes() == (c2 ** 0.67).tobytes()


def test_the_factor_refuses_other_genes_bad_values_an_undefined_c_and_other_points():
    m = np.array([2.0, 0.5, 1.0])
    with pytest.raises(bk.Refused, match="same genes"):
        bk.c2_gonly_factor(m, np.array([2.0, 0.5]), INF)
    with pytest.raises(bk.Refused, match="finite and positive"):
        bk.c2_gonly_factor(m, np.array([2.0, 0.0, 1.0]), INF)
    with pytest.raises(bk.Refused, match="c_inf is undefined"):
        bk.c2_gonly_factor(m, np.ones(3), INF)
    for g in (2.0, 0.0, -INF, float("nan")):
        with pytest.raises(bk.Refused, match="not this line's"):
            bk.c2_gonly_factor(m, m[::-1], g)
    for s in (0.0, -0.67, INF, float("nan")):
        with pytest.raises(bk.Refused, match="size scale"):
            bk.c2_gonly_factor(m, m[::-1], INF, s)


def test_combine_and_the_c2_gata4_weight_are_unchanged():
    phi, rec = bk.combine("c1", 1.0, PHI_M, PHI_G)
    assert np.array_equal(phi, th.C1_HAND) and rec["gata4_weight"] == 1.0
    lm, lg = np.log(PHI_M), np.log(PHI_G)
    for x in (1.0, 1.5, 2.0):
        phi, rec = bk.combine("c2", x, PHI_M, PHI_G)
        c = np.std(x * lm) / np.std(x * lm + lg)
        assert rec["c"] == pytest.approx(c, rel=1e-13) and np.array_equal(phi, (PHI_M ** (x * c)) * (PHI_G ** c))
    with pytest.raises(bk.Refused, match="not finite"):
        bk.c2_gata4_weighted(PHI_M, PHI_G, INF)  # line 3's own function still takes no infinite weight


# ------------------------------------------------------------------------------------------------------- the checks


@pytest.mark.parametrize("key", list(TABLES))
@pytest.mark.parametrize("name", list(POINTS))
def test_the_range_guard_by_hand(key, name):
    m, t = TABLES[key]
    g, s = POINTS[name]
    r = bk.gonly_guards(m, t, g, s, np.ones(G), lift={"row": "R44"})  # the den check lifted: the range alone gates
    lp = ln_factor_by_hand(m, t, g, s)
    c2 = np.exp(ln_factor_by_hand(m, t, 1.0))
    rg = r["range"]
    assert rg["factor_min"] == pytest.approx(float(np.exp(lp).min()), rel=1e-12)
    assert rg["factor_max"] == pytest.approx(float(np.exp(lp).max()), rel=1e-12)
    assert rg["c2_min"] == pytest.approx(float(c2.min()), rel=1e-12)
    assert rg["c2_max"] == pytest.approx(float(c2.max()), rel=1e-12)
    assert rg["bound_min"] == rg["c2_min"] / 2.0 and rg["bound_max"] == rg["c2_max"] * 2.0
    assert rg["passed"] is (rg["factor_min"] >= rg["bound_min"] and rg["factor_max"] <= rg["bound_max"])
    assert rg["passed"] is RANGE_OK[key][name] and r["passed"] is RANGE_OK[key][name]
    assert set(r) == {"g", "size_scale", "c", "powers_applied", "range", "den_check", "passed", "factor_sha256"}
    phi, _ = bk.c2_gonly_factor(m, t, g, s)
    assert r["factor_sha256"] == hashlib.sha256(np.ascontiguousarray(phi).tobytes()).hexdigest()


@pytest.mark.parametrize("edge,ok", [("max_at_bound", True), ("max_beyond", False), ("min_at_bound", True),
                                     ("min_beyond", False)])
def test_the_range_bound_is_inclusive_at_twice_and_half(edge, ok, monkeypatch):
    """The applied factor at the bounds exactly (the factor stubbed): max = 2 x C2's max passes, one ulp beyond fails;
    min = C2's min / 2 passes, one ulp below fails."""
    c2, _ = bk.combine("c2", 1.0, *ORTH)
    f = np.ones(G)
    if edge.startswith("max"):
        f[0] = 2.0 * c2.max() if edge == "max_at_bound" else np.nextafter(2.0 * c2.max(), INF)
    else:
        f[0] = c2.min() / 2.0 if edge == "min_at_bound" else np.nextafter(c2.min() / 2.0, 0.0)
    monkeypatch.setattr(bk, "c2_gonly_factor", lambda *a, **k: (f, {"g": "inf", "size_scale": 0.67, "c": 1.0,
                                                                     "mab21l2_exponent": 0.0, "gata4_weight": 0.67}))
    r = bk.gonly_guards(*ORTH, INF, 0.67, np.ones(G), lift={"row": "R44"})
    assert r["range"]["passed"] is ok and r["passed"] is ok


@pytest.mark.parametrize("name", list(POINTS))
def test_the_den_check_by_hand_is_line_3_s(name):
    """Line 3's check (c2g_guard's figures): the file's den-weighted total against C2's factor's, within 5%; here by
    hand on one-hot and spread dens over tables where it passes (MILD) and fails (SHIFT)."""
    g, s = POINTS[name]
    for key, dens in (("mild", [np.ones(G), np.linspace(0.1, 3.0, G)]), ("shift", [np.ones(G),
                                                                                  np.linspace(3.0, 0.1, G)])):
        m, t = TABLES[key]
        phi, _ = bk.c2_gonly_factor(m, t, g, s)
        c2, _ = bk.combine("c2", 1.0, m, t)
        for den in dens:
            r = bk.gonly_guards(m, t, g, s, den)
            d = r["den_check"]
            tot = den.sum()
            assert d["den_total_phi_m"] == pytest.approx(den @ m / tot, rel=1e-14)
            assert d["den_total_c2"] == pytest.approx(den @ c2 / tot, rel=1e-14)
            assert d["den_total_factor"] == pytest.approx(den @ phi / tot, rel=1e-14)
            assert d["den_shift_vs_c2"] == pytest.approx((den @ phi) / (den @ c2) - 1.0, rel=1e-12, abs=1e-15)
            assert d["den_shift_c2_vs_phi_m"] == pytest.approx((den @ c2) / (den @ m) - 1.0, rel=1e-12, abs=1e-15)
            assert d["den_shift_factor_vs_phi_m"] == pytest.approx((den @ phi) / (den @ m) - 1.0, rel=1e-12,
                                                                    abs=1e-15)
            assert d["max_den_shift"] == 0.05 and d["vs_phi_m"] == "recorded; gates nothing"
            assert d["within"] is (abs(d["den_shift_vs_c2"]) <= 0.05) is (key == "mild")
            assert d["passed"] is d["within"] and d["lifted_by"] is None and d["limit"]["status"] == "gates"
            assert r["passed"] is (key == "mild") and r["range"]["passed"]
            assert d["rule"] == bk.GONLY_DEN_RULE
    with pytest.raises(bk.Refused, match="den"):
        bk.gonly_guards(*MILD, g, s, np.ones(G - 1))
    with pytest.raises(bk.Refused, match="den"):
        bk.gonly_guards(*MILD, g, s, -np.ones(G))


@pytest.mark.parametrize("shift,ok", [(0.0499, True), (-0.0499, True), (0.0, True), (0.0501, False),
                                      (-0.0501, False), (0.5, False)])
def test_the_den_bound_is_5_percent_either_way_and_a_lift_makes_it_gate_nothing(shift, ok, monkeypatch):
    """The bound on the shift against C2's total, both signs, at the edges (C2's factor stubbed as the file's over
    1 + shift); with a lifting record the check records the same figures and passes."""
    real_factor = bk.c2_gonly_factor

    def scaled(kind, x, pm, pg):
        file_phi, _ = real_factor(pm, pg, INF, 1.0)
        return file_phi / (1.0 + shift), {"c": 1.0}

    monkeypatch.setattr(bk, "combine", scaled)
    r = bk.gonly_guards(*ORTH, INF, 1.0, np.ones(G))
    d = r["den_check"]
    assert d["den_shift_vs_c2"] == pytest.approx(shift, abs=1e-12) and d["within"] is ok and d["passed"] is ok
    assert r["passed"] is ok and r["range"]["passed"]
    lift = {"row": "R44", "option": f"`{OPT['F1']}`"}
    lr = bk.gonly_guards(*ORTH, INF, 1.0, np.ones(G), lift=lift)["den_check"]
    assert lr["passed"] is True and lr["within"] is ok and lr["lifted_by"] == lift
    assert lr["den_shift_vs_c2"] == d["den_shift_vs_c2"]
    assert lr["limit"] == {"max_abs_den_shift": 0.05, "status": "lifted by row R44: the 5% check gates nothing for it"}


def test_a_failed_check_says_which():
    r = bk.gonly_guards(*WIDE, INF, 1.0, np.ones(G), lift={"row": "R44"})
    rg = r["range"]
    assert bk.gonly_failure(r) == (f"the range guard: the factor spans [0.0625, 16], not within "
                                   f"[{rg['bound_min']:.4g}, {rg['bound_max']:.4g}] (half C2's min, twice C2's max): "
                                   "not built")
    assert rg["bound_max"] == pytest.approx(2.0 * 16.0 ** 2.0 ** -0.5, rel=1e-12)
    d = bk.gonly_guards(*SHIFT, INF, 1.0, np.ones(G))
    says = bk.gonly_failure(d)
    assert says.startswith("the 5% den check: the file's den-weighted total is +") and says.endswith(
        "against C2's factor's, beyond 5%, and no --den-ruling row lifts it for this file: not built")
    both = bk.gonly_guards(*SPIKE, INF, 1.0, np.ones(G))
    if not both["den_check"]["passed"]:
        assert "the range guard" in bk.gonly_failure(both) and "; the 5% den check" in bk.gonly_failure(both)


# ------------------------------------------------------------------------------------------------- the ruling rows


@pytest.mark.parametrize("rows,ok", [
    (ROW, {"F1", "F2", "F3"}),
    (ROW_NEW, {"F1", "F2", "F3"}),
    (row_only("F1"), {"F1"}),
    (row_only("F2"), {"F2"}),
    (row_only("F3"), {"F3"}),  # F3's entry opens F3 only, not F1, whose string it contains
    (row_only("F1", "F2"), {"F1", "F2"}),
    (ROW_OLD, set()),  # the old lines' entries open none of the new files
    (ROW + ROW, set()),  # two rows R22
    ("", set()),  # no row R22
    (ROW.replace("| R22 (the captain's order", "| R22 (the captain's ruling"), set()),
    (ROW.replace("ruling, no run", "a run"), set()),
    (ROW.replace("| R22 (", "| R22b ("), set()),
    (ROW.replace("`--captain-ruling R22`", "R22"), set()),
    (ROW.replace(f"`{OPT['F1']}`", OPT["F1"]), {"F2", "F3"}),  # F1's without its backticks
    (ROW.replace(f"`{OPT['F1']}`", "`--gata4-exponent inf`"), {"F2", "F3"}),  # without its form
    (ROW.replace(f"`{OPT['F1']}`", f"`a`{OPT['F1']}`b`"), {"F2", "F3"}),  # in the text, not as one whole string
    (ROW.replace(f"`{OPT['F3']}`", "`--gata4-exponent inf --size-scale 0.67`"), {"F1", "F2"}),
    (ROW.replace(f"`{OPT['F2']}`", "`--size-scale 0.67`"), {"F1", "F3"}),
    (ROW.replace(f"`{OPT['F2']}`", "`--form c2-s2 --size-scale 0.670`"), {"F1", "F3"}),
    (ROW_OLD.replace(") | T3 | x |", ") | T3 | " + " ".join(f"`{o}`" for o in OPT.values()) + " |"), set()),
])
def test_the_ruling_row_check(tmp_path, rows, ok):
    p = ledger(tmp_path, rows)
    for name in POINTS:
        if name in ok:
            rec = bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, **OPTS[name])
            assert rec == {"row": "R22", "file": FILES[name], "name": name, "option": f"`{OPT[name]}`",
                           "g": "inf" if POINTS[name][0] == INF else 1.0, "size_scale": POINTS[name][1],
                           "den_ruling": None, "ledger": "experiments/harness/TRIES.md",
                           "ledger_sha256": t3d.sha256(p)}
            assert bk.gonly_point(rec) == POINTS[name]
        else:
            with pytest.raises(bk.Refused, match="no single row R22"):
                bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, **OPTS[name])


def test_the_new_entries_do_not_open_the_old_lines_and_the_old_entries_open_them_still(tmp_path):
    new_only = ledger(tmp_path, ROW_NEW, name="new.md")
    with pytest.raises(bk.Refused, match="no single row R22"):
        bk.check_gata4("c1-s2", 2.0, "R22", new_only)
    with pytest.raises(bk.Refused, match="no single row R22"):
        bk.check_c2_gata4("c2-s2", 2.0, "R22", new_only)
    both = ledger(tmp_path, ROW, name="both.md")
    assert bk.check_gata4("c1-s2", 2.0, "R22", both)["option"] == "`--gata4-exponent 2.0`"
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", both)["option"] == "`--form c2-s2 --gata4-exponent 2.0`"


def test_the_row_check_is_build_outside_form_s_plus_the_whole_string(tmp_path, monkeypatch):
    bo = bk._outside_form()
    seen = []
    real = bo.captain_order_names
    monkeypatch.setattr(bo, "captain_order_names", lambda row, option, path=None: seen.append((row, option, path))
                        or real(row, option, path))
    p = ledger(tmp_path)
    for name in POINTS:
        seen.clear()
        assert bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, **OPTS[name])["name"] == name
        assert seen == [("R22", f"`{OPT[name]}`", p)]
    # captain_order_names alone would pass the odd cell; the whole-string match refuses it
    odd = ledger(tmp_path, ROW.replace(f"`{OPT['F1']}`", f"`a`{OPT['F1']}`b`"), name="odd.md")
    assert real("R22", f"`{OPT['F1']}`", odd)
    with pytest.raises(bk.Refused, match="no single row R22"):
        bk.check_c2_gonly("c2-s2", INF, None, "R22", odd)
    assert bk.check_c2_gonly("c2-s2", INF, 0.67, "R22", odd)["name"] == "F3"  # the later strings are unaffected


@pytest.mark.parametrize("name", list(POINTS))
def test_a_lifting_row_that_names_the_file_lifts_it(tmp_path, name):
    p = ledger(tmp_path, ROW + lift_row(name))
    rec = bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, den_ruling="R44", **OPTS[name])
    assert rec["den_ruling"] == {"row": "R44", "option": f"`{OPT[name]}`", "words": "the 5% check gates nothing for it",
                                 "ledger": "experiments/harness/TRIES.md", "ledger_sha256": t3d.sha256(p)}
    assert bk.check_den_lift(name, "R44", p) == rec["den_ruling"]


@pytest.mark.parametrize("lift,lifted", [
    (lift_row("F1", "F3"), {"F1", "F3"}),
    (lift_row("F2"), {"F2"}),
    (lift_row("F3"), {"F3"}),  # F3's string does not lift F1, whose string it contains
    (lift_row("F1", "F2", "F3"), {"F1", "F2", "F3"}),
    (lift_row("F1", "F2", "F3", words="the 5% check gates nothing for them"), set()),
    (lift_row("F1", "F2", "F3", words="the 5% check gates nothing"), set()),
    (lift_row("F1", "F2", "F3", ruling="a note"), set()),  # not a ruling
    (lift_row("F1", "F2", "F3") * 2, set()),  # two rows R44
    (lift_row("F1", "F2", "F3", row="R45"), set()),  # no row R44
    (lift_row("F1", "F2", "F3").replace("`", ""), set()),  # the strings without their backticks
    ("| R44 (the captain's yes of a day; ruling, no run; for " + " and ".join(f"`{o}`" for o in OPT.values())
     + ") | T3 | the 5% check gates nothing for it |\n", set()),  # the words outside the first cell
    ("| R44 (the captain's yes of a day; ruling, no run; the 5% check gates nothing for it) | T3 | "
     + " ".join(f"`{o}`" for o in OPT.values()) + " |\n", set()),  # the strings outside the first cell
])
def test_the_lifting_row_check(tmp_path, lift, lifted):
    p = ledger(tmp_path, ROW + lift)
    for name in POINTS:
        if name in lifted:
            assert bk.check_den_lift(name, "R44", p)["row"] == "R44"
        else:
            with pytest.raises(bk.Refused, match="--den-ruling R44 is refused: .*no single row R44"):
                bk.check_den_lift(name, "R44", p)
            with pytest.raises(bk.Refused, match="no single row R44"):
                bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, den_ruling="R44", **OPTS[name])


@pytest.mark.parametrize("value", ["r44", "R", "R44a", "44", "", " R44", "R44 ", "R-1", "R4 4", None, 44, True])
def test_a_lifting_row_id_other_than_r_and_digits_is_refused(tmp_path, value):
    p = ledger(tmp_path, ROW + lift_row("F1", "F2", "F3"))
    with pytest.raises(bk.Refused, match="takes a ruling row's id"):
        bk.check_den_lift("F1", value, p)


# ------------------------------------------------------------------------------------------------------ refusals


def test_without_the_options_nothing_of_the_line_is_read(monkeypatch):
    monkeypatch.setattr(bk, "_outside_form", lambda: pytest.fail("read on a path without the line"))
    for form in [*bk.FORMS, "c1-s1x1.5", "", None]:
        for kw in ({}, {"captain_ruling": "R22"}, {"captain_ruling": "R36"}, {"gata4_exponent": 2.0},
                   {"gata4_exponent": 3.0, "captain_ruling": "R22"}, {"composition_from": "allowed-s1"},
                   {"gata4_exponent": -INF}, {"gata4_exponent": float("nan")}):
            assert bk.check_c2_gonly(form, **kw) is None and not bk.is_gonly_command(form, kw.get("gata4_exponent"))
    for form in [f for f in bk.FORMS if f != "c2-s2"]:
        assert bk.check_c2_gonly(form, INF, None, "R22") is None  # inf elsewhere: refused as before (only 2.0)
    assert bk.is_gonly_command("c2-s2", INF) and bk.is_gonly_command("c2-s2", "inf")
    assert bk.is_gonly_command("c1-s2", None, 0.67) and not bk.is_gonly_command("c2-s2", True)


@pytest.mark.parametrize("form,kw", [("c2-s2", {}), ("c2-s2", {"captain_ruling": "R22"}),
                                     ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}),
                                     ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}),
                                     ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36"}),
                                     ("identity-s2", {}), ("c1-s2", {"gata4_exponent": INF, "captain_ruling": "R22"})])
def test_the_lifting_option_on_any_other_command_is_refused(tmp_path, form, kw):
    with pytest.raises(bk.Refused, match="--den-ruling is taken only with the options of F1, F2 or F3"):
        bk.check_c2_gonly(form, den_ruling="R44", tries=ledger(tmp_path, ROW + lift_row("F1", "F2", "F3")), **kw)


@pytest.mark.parametrize("value", [0.5, 1.0, 1.5, 2.0, 0.66, 0.68, 0.6700001, 0.0, -0.67, float("nan"), INF, "x", "",
                                   True, False, "0,67"])
@pytest.mark.parametrize("g", [None, INF])
def test_any_other_size_scale_is_refused(tmp_path, value, g):
    with pytest.raises(bk.Refused, match="only 0.67"):
        bk.check_c2_gonly("c2-s2", g, value, "R22", ledger(tmp_path))


@pytest.mark.parametrize("value", [2.0, 3.0, 1.0, 0.0, -INF, float("nan"), "x", True])
def test_the_size_scale_beside_another_weight_is_refused(tmp_path, value):
    with pytest.raises(bk.Refused, match="taken alone .* or with --gata4-exponent inf"):
        bk.check_c2_gonly("c2-s2", value, 0.67, "R22", ledger(tmp_path))


@pytest.mark.parametrize("form", [f for f in bk.FORMS if f != "c2-s2"])
def test_the_size_scale_on_another_form_is_refused(tmp_path, form):
    for g in (None, INF, 2.0):
        with pytest.raises(bk.Refused, match="--size-scale is taken on --form c2-s2 only"):
            bk.check_c2_gonly(form, g, 0.67, "R22", ledger(tmp_path))


@pytest.mark.parametrize("ruling", ["R21", "R30", "R41", "r22", "R22 ", " R22", "R22b", "", "R36"])
@pytest.mark.parametrize("name", list(POINTS))
def test_any_other_ruling_is_refused(tmp_path, ruling, name):
    with pytest.raises(bk.Refused, match="only row R22"):
        bk.check_c2_gonly("c2-s2", captain_ruling=ruling, tries=ledger(tmp_path), **OPTS[name])


@pytest.mark.parametrize("name", list(POINTS))
def test_either_option_without_the_ruling_or_beside_the_composition_option_is_refused(tmp_path, name):
    with pytest.raises(bk.Refused, match="only with --captain-ruling R22"):
        bk.check_c2_gonly("c2-s2", tries=ledger(tmp_path), **OPTS[name])
    with pytest.raises(bk.Refused, match="takes no --composition-from"):
        bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=ledger(tmp_path), composition_from="allowed-s1",
                          **OPTS[name])


@pytest.mark.parametrize("form,kw,says", [
    ("c2-s2", {"size_scale": 0.5, "captain_ruling": "R22"}, "only 0.67"),
    ("c2-s2", {"size_scale": 0.67, "gata4_exponent": 2.0, "captain_ruling": "R22"}, "taken alone"),
    ("c2-s2", {"gata4_exponent": INF}, "only with --captain-ruling R22"),
    ("c2-s2", {"size_scale": 0.67}, "only with --captain-ruling R22"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R30"}, "only row R22"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R36"}, "only row R22"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R36", "composition_from": "allowed-s1"}, "only row R22"),
    ("c2-s2", {"size_scale": 0.67, "captain_ruling": "R36", "composition_from": "allowed-s1"}, "only row R22"),
    ("c2-s2", {"size_scale": 0.67, "captain_ruling": "R22", "composition_from": "allowed-s1"}, "no --composition"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R22", "rows": ROW_OLD}, "no single row R22"),
    ("c2-s2", {"size_scale": 0.67, "captain_ruling": "R22", "rows": row_only("F1")}, "no single row R22"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R22", "rows": row_only("F3")}, "no single row R22"),
    ("c2-s2", {"gata4_exponent": INF, "size_scale": 0.67, "captain_ruling": "R22", "rows": row_only("F1", "F2")},
     "no single row R22"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R22", "den_ruling": "R44", "rows": ROW + lift_row("F3")},
     "no single row R44"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R22", "den_ruling": "R4x"}, "takes a ruling row's id"),
    ("c2-s2", {"captain_ruling": "R22", "den_ruling": "R44"}, "--den-ruling is taken only"),
    ("c1-s2", {"size_scale": 0.67, "captain_ruling": "R22"}, "on --form c2-s2 only"),
    ("identity-s2", {"size_scale": 0.67, "captain_ruling": "R22"}, "on --form c2-s2 only"),
    ("c1-s2", {"gata4_exponent": INF, "captain_ruling": "R22"}, "only 2.0"),  # inf elsewhere: as before
    ("c2-s1", {"gata4_exponent": INF, "captain_ruling": "R22"}, "only 2.0"),
    ("c2-s2", {"gata4_exponent": -INF, "captain_ruling": "R22"}, "only 2.0"),  # -inf: line 3's refusal, as before
])
def test_the_builder_refuses_and_builds_nothing(world, tmp_path, monkeypatch, form, kw, says):
    for name in (*LINE_FUNCS[:-1], "run_c2g_guard", "c2g_loader", "run_gata4_figures"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on a refusal"))
    kw = dict(kw)
    t = ledger(tmp_path, kw.pop("rows", ROW))
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match=says):
        build(world, out, form, tries=t, **kw)
    assert nothing_written(out) and loader_restored()


# ------------------------------------------------------------------------------------------- the files at the build


@pytest.mark.parametrize("name", list(POINTS))
def test_each_file_is_c2_s_command_with_its_factor_and_records_both_checks(world, runs, tmp_path, monkeypatch, name):
    r = use(runs, monkeypatch, "mild")
    spy = Spy(monkeypatch)
    c2 = build(world, tmp_path / "c2.h5ad", run_dir=r.run)
    out = tmp_path / FILES[name]
    side = build_f(world, tmp_path, name, out, run_dir=r.run)
    assert len(spy.seen) == 2
    want, rec = bk.c2_gonly_factor(r.m, r.g, *POINTS[name])
    assert np.array_equal(spy.seen[1], want) and np.array_equal(spy.seen[0], bk.combine("c2", 1.0, r.m, r.g)[0])
    # C2's command with its own options: spec, inputs, streams, composition, seed, count and board are C2's
    assert side["spec"] == c2["spec"] and side["inputs"] == c2["inputs"] and side["spec_source"] == c2["spec_source"]
    assert side["inputs"]["rng_streams"] == ["submission:barred:s2", "submission:barred:s2:assemble"]
    assert side["construction"]["composition"] == c2["construction"]["composition"]
    assert side["construction"]["factor_exponent"] == 1.0 and side["output"]["check"]["ok"]
    cb, cb2 = side["combined_build"], c2["combined_build"]
    assert {k: v for k, v in cb.items() if k not in ("combination", "c2_gonly_line")} == {
        k: v for k, v in cb2.items() if k != "combination"}
    assert cb["form"] == "c2-s2" and cb["combination"]["c"] == rec["c"]
    line = cb["c2_gonly_line"]
    assert line["captain_ruling"]["row"] == "R22" and line["captain_ruling"]["option"] == f"`{OPT[name]}`"
    assert line["file"] == FILES[name] and line["g"] == ("inf" if POINTS[name][0] == INF else 1.0)
    assert line["size_scale"] == POINTS[name][1] and line["c"] == rec["c"] and line["line"] == bk.GONLY_LINE
    assert line["powers_applied"] == {bk.TABLE_M: rec["mab21l2_exponent"], bk.TABLE_G: rec["gata4_weight"]}
    assert line["den_ruling"] is None and line["captain_ruling"]["den_ruling"] is None
    sha = hashlib.sha256(np.ascontiguousarray(want).tobytes()).hexdigest()
    assert line["factor_sha256_applied"] == line["guards"]["factor_sha256"] == sha
    gd = line["guards"]
    assert "genes" not in gd and gd["n_genes"] == G and gd["tables_sha256"] == r.pins
    assert gd["passed"] and gd["range"]["passed"] and gd["den_check"]["passed"] and gd["den_check"]["within"]
    assert gd["den_check"]["limit"]["status"] == "gates" and "distinctness" not in gd
    want_g = bk.gonly_guards(r.m, r.g, *POINTS[name], den_of(world))

    def as_json(o):
        return json.loads(json.dumps(t3d.json_safe(o)))

    assert as_json({k: gd[k] for k in want_g}) == as_json(want_g)
    assert np.array_equal(np.array(gd["den"]), den_of(world)) and gd["wt_sha256"] == t3d.sha256(world.wt)
    on_disk = json.loads(out.with_suffix(".json").read_text())["combined_build"]["c2_gonly_line"]
    assert on_disk["file"] == FILES[name] and on_disk["guards"]["range"] == gd["range"]
    assert json.loads(Path(f"{out}.stamp.json").read_text())["clean"] is True
    assert not any("cell" in k.lower() for k in line) and loader_restored()
    assert "cos" not in json.dumps(on_disk).lower().replace("cost", "")


@pytest.mark.parametrize("name", list(POINTS))
def test_each_file_draws_c2_s_cells_differs_in_expression_and_builds_one_sha256(world, runs, tmp_path, monkeypatch,
                                                                                 name):
    import anndata as ad

    r = use(runs, monkeypatch, "mild")
    build(world, tmp_path / "c2.h5ad", run_dir=r.run)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    sa = build_f(world, tmp_path, name, tmp_path / "a" / FILES[name], run_dir=r.run)
    sb = build_f(world, tmp_path, name, tmp_path / "b" / FILES[name], run_dir=r.run)
    assert (tmp_path / "a" / FILES[name]).read_bytes() == (tmp_path / "b" / FILES[name]).read_bytes()
    assert sa["output"]["sha256"] == sb["output"]["sha256"]
    a, c = ad.read_h5ad(tmp_path / "a" / FILES[name]), ad.read_h5ad(tmp_path / "c2.h5ad")
    assert list(a.obs_names) == list(c.obs_names) and a.obs.equals(c.obs)
    assert np.array_equal(a.obsm["spatial_3D"], c.obsm["spatial_3D"])
    assert not np.array_equal(np.asarray(a.X), np.asarray(c.X))


@pytest.mark.parametrize("key,name,says", [("wide", "F1", "the range guard"), ("spike", "F1", "the range guard"),
                                           ("spike", "F3", "the range guard"), ("shift", "F1", "the 5% den check"),
                                           ("shift", "F2", "the 5% den check"), ("shift", "F3", "the 5% den check")])
def test_a_failed_check_stops_the_file_before_any_cell_and_writes_nothing(world, runs, tmp_path, monkeypatch, key,
                                                                          name, says):
    r = use(runs, monkeypatch, key)
    monkeypatch.setattr(bk, "gonly_loader", lambda *a, **k: pytest.fail("the build was reached past a failed check"))
    monkeypatch.setattr(bk.bs, "construct", lambda *a, **k: pytest.fail("a cell was drawn past a failed check"))
    out = tmp_path / FILES[name]
    with pytest.raises(bk.GonlyGuardFailed, match=says):
        build_f(world, tmp_path, name, out, run_dir=r.run)
    assert nothing_written(out) and loader_restored()


@pytest.mark.parametrize("name", list(POINTS))
def test_the_lifting_row_lets_a_file_the_5_percent_check_stops_build(world, runs, tmp_path, monkeypatch, name):
    """SHIFT: each file's den total is far from C2's, so the check stops it; with --den-ruling naming a row that names
    the file and says the words, it builds and records the figures, which gate nothing."""
    r = use(runs, monkeypatch, "shift")
    out = tmp_path / FILES[name]
    side = build_f(world, tmp_path, name, out, rows=ROW + lift_row(name), run_dir=r.run, den_ruling="R44")
    line = side["combined_build"]["c2_gonly_line"]
    d = line["guards"]["den_check"]
    assert d["passed"] and not d["within"] and abs(d["den_shift_vs_c2"]) > 0.05
    assert d["limit"]["status"] == "lifted by row R44: the 5% check gates nothing for it"
    assert d["lifted_by"]["row"] == "R44" and d["lifted_by"]["option"] == f"`{OPT[name]}`"
    assert line["den_ruling"] == d["lifted_by"] == line["captain_ruling"]["den_ruling"]
    assert line["guards"]["range"]["passed"] and side["output"]["check"]["ok"]
    # a row that lifts another file only is refused for this one: exit 2, nothing written
    other = {"F1": "F2", "F2": "F3", "F3": "F1"}[name]
    out2 = tmp_path / f"other_{FILES[name]}"
    with pytest.raises(bk.Refused, match="no single row R44"):
        build_f(world, tmp_path, name, out2, rows=ROW + lift_row(other), run_dir=r.run, den_ruling="R44")
    assert nothing_written(out2)


def test_the_lifting_row_does_not_lift_the_range_guard(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "wide")
    out = tmp_path / FILES["F1"]
    with pytest.raises(bk.GonlyGuardFailed, match="the range guard"):
        build_f(world, tmp_path, "F1", out, rows=ROW + lift_row("F1"), run_dir=r.run, den_ruling="R44")
    assert nothing_written(out)


def test_a_factor_other_than_the_checks_exits_5_and_writes_nothing(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    real = bk.run_gonly_guards
    monkeypatch.setattr(bk, "run_gonly_guards", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.FactorMismatch, match="not the one the guards describe"):
        build_f(world, tmp_path, "F1", out, run_dir=r.run)
    assert nothing_written(out)
    monkeypatch.setattr(bk, "run_gonly_guards", lambda *a, **k: real(*a, **k) | {"genes": ["Z"]})
    with pytest.raises(bk.FactorMismatch, match="genes"):
        build_f(world, tmp_path, "F2", out, run_dir=r.run)
    assert nothing_written(out) and loader_restored()


def test_the_checks_read_the_pinned_tables_only(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    monkeypatch.setattr(bk, "TABLE_SHA256", {**r.pins, bk.TABLE_G: "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match="pinned table"):
        build_f(world, tmp_path, "F3", out, run_dir=r.run)
    assert nothing_written(out)


def test_an_existing_output_is_refused_before_the_checks(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    out = tmp_path / FILES["F1"]
    out.write_bytes(b"x")
    monkeypatch.setattr(bk, "run_gonly_guards", lambda *a, **k: pytest.fail("the checks ran past an existing file"))
    with pytest.raises(bk.Refused, match="exists; pass --overwrite"):
        build_f(world, tmp_path, "F1", out, run_dir=r.run)
    assert out.read_bytes() == b"x"


# ------------------------------------------------------------------------------- every existing path as it was


OLD_SETS = [(f, {}) for f in bk.FORMS] + [("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"})]


@pytest.mark.parametrize("form,kw", OLD_SETS, ids=[f"{f}{'-g2' if kw else ''}" for f, kw in OLD_SETS])
def test_every_existing_command_never_touches_the_line(world, tmp_path, monkeypatch, form, kw):
    t = ledger(tmp_path)
    extra = dict(kw, tries=t) if kw else {}
    plain = build(world, tmp_path / "plain.h5ad", form, **extra)
    for name in LINE_FUNCS:
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on {form}"))
    trapped = build(world, tmp_path / "trapped.h5ad", form, size_scale=None, den_ruling=None, **extra)
    assert (tmp_path / "plain.h5ad").read_bytes() == (tmp_path / "trapped.h5ad").read_bytes()
    assert plain["combined_build"] == trapped["combined_build"] and "c2_gonly_line" not in plain["combined_build"]


def test_line_3_s_file_never_touches_the_line(world, tmp_path, monkeypatch):
    """`--form c2-s2 --gata4-exponent 2.0` on a run that passes line 3's own 5% guard (the quarter table)."""
    run = world.d / "run_quarter_gonly"
    if not run.exists():
        th.write_run(run, world.wt, {bk.TABLE_M: PHI_M, bk.TABLE_G: PHI_G ** 0.25})
    monkeypatch.setattr(bk, "TABLE_SHA256", {k: t3d.sha256(run / "factors" / f"{k}.csv")
                                             for k in (bk.TABLE_M, bk.TABLE_G)})
    t = ledger(tmp_path)
    plain = build(world, tmp_path / "g2.h5ad", run_dir=run, gata4_exponent=2.0, captain_ruling="R22", tries=t)
    for name in LINE_FUNCS:
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on line 3"))
    trapped = build(world, tmp_path / "g2b.h5ad", run_dir=run, gata4_exponent=2.0, captain_ruling="R22", tries=t)
    assert (tmp_path / "g2.h5ad").read_bytes() == (tmp_path / "g2b.h5ad").read_bytes()
    assert trapped["combined_build"]["c2_gata4_line"]["gata4_weight_g"] == 2.0
    assert plain["combined_build"] == trapped["combined_build"] and "c2_gonly_line" not in trapped["combined_build"]


# --------------------------------------------------------------------------------------------------------- CLI


def test_the_cli(world, runs, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path, ROW + lift_row("F1", "F2", "F3")))
    boards = th.cli_boards(tmp_path)
    r = use(runs, monkeypatch, "mild")

    def args(run, *extra, out):
        return ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
                str(boards), "--run-dir", str(run), "--seed", str(SEED), "--threads", "2", *extra, "--out", str(out)]

    runner = CliRunner()
    for extra in (["--form", "c2-s2", "--size-scale", "0.5", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--size-scale", "0.67"],
                  ["--form", "c2-s2", "--gata4-exponent", "inf"],
                  ["--form", "c2-s2", "--gata4-exponent", "inf", "--captain-ruling", "R41"],
                  ["--form", "c2-s2", "--gata4-exponent", "2.0", "--size-scale", "0.67", "--captain-ruling", "R22"],
                  ["--form", "c1-s2", "--size-scale", "0.67", "--captain-ruling", "R22"],
                  ["--form", "c1-s2", "--gata4-exponent", "inf", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--gata4-exponent=-inf", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--size-scale", "0.67", "--captain-ruling", "R36", "--composition-from",
                   "allowed-s1"],
                  ["--form", "c2-s2", "--size-scale", "0.67", "--captain-ruling", "R22", "--composition-from",
                   "allowed-s1"],
                  ["--form", "c2-s2", "--den-ruling", "R44"],
                  ["--form", "c2-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R22", "--den-ruling", "R44"],
                  ["--form", "c2-s2", "--gata4-exponent", "inf", "--captain-ruling", "R22", "--den-ruling", "R4a"],
                  ["--form", "c2-s2", "--gata4-exponent", "inf", "--captain-ruling", "R22", "--den-ruling", "R45"]):
        out = tmp_path / "r.h5ad"
        res = runner.invoke(bk.app, args(r.run, *extra, out=out))
        assert res.exit_code == 2 and "refused" in res.output, (extra, res.output)
        assert nothing_written(out)
    lines = {"F1": ["--gata4-exponent", "inf"], "F2": ["--size-scale", "0.67"],
             "F3": ["--gata4-exponent", "inf", "--size-scale", "0.67"]}
    for name, opts in lines.items():
        ok = tmp_path / FILES[name]
        res = runner.invoke(bk.app, args(r.run, "--form", "c2-s2", *opts, "--captain-ruling", "R22", out=ok))
        assert res.exit_code == 0, res.output
        assert (f"row R22's options on C2's command: {FILES[name]}" in res.output
                and "the range guard passed; the 5% den check passed" in res.output)
        line = json.loads(ok.with_suffix(".json").read_text())["combined_build"]["c2_gonly_line"]
        assert line["file"] == FILES[name] and line["captain_ruling"]["name"] == name
    # a check fails: exit 6, nothing written (WIDE: F1's range; SHIFT: the 5% den check)
    for key, says in (("wide", "the range guard"), ("shift", "the 5% den check")):
        use(runs, monkeypatch, key)
        bad = tmp_path / f"stop_{key}.h5ad"
        res = runner.invoke(bk.app, args(runs[key].run, "--form", "c2-s2", *lines["F1"], "--captain-ruling", "R22",
                                         out=bad))
        assert res.exit_code == 6 and "refused (the range guard or the 5% den check" in res.output, res.output
        assert says in res.output and nothing_written(bad)
    # the lifting row: SHIFT's F2 builds with --den-ruling R44
    use(runs, monkeypatch, "shift")
    lifted = tmp_path / f"lifted_{FILES['F2']}"
    res = runner.invoke(bk.app, args(runs["shift"].run, "--form", "c2-s2", *lines["F2"], "--captain-ruling", "R22",
                                     "--den-ruling", "R44", out=lifted))
    assert res.exit_code == 0 and "the 5% den check gates nothing for it (row R44)" in res.output, res.output
    # the factor is not the checks': exit 5, nothing written
    use(runs, monkeypatch, "mild")
    real = bk.run_gonly_guards
    monkeypatch.setattr(bk, "run_gonly_guards", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    mis = tmp_path / "mismatch.h5ad"
    res = runner.invoke(bk.app, args(r.run, "--form", "c2-s2", *lines["F3"], "--captain-ruling", "R22", out=mis))
    assert res.exit_code == 5 and nothing_written(mis)
    monkeypatch.setattr(bk, "run_gonly_guards", real)
    # the ruling row naming the old lines only: refused
    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path, ROW_OLD, name="old.md"))
    for opts in lines.values():
        res = runner.invoke(bk.app, args(r.run, "--form", "c2-s2", *opts, "--captain-ruling", "R22",
                                         out=tmp_path / "n.h5ad"))
        assert res.exit_code == 2 and "no single row R22" in res.output and nothing_written(tmp_path / "n.h5ad")
