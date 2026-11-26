"""Row R76's five files in experiments/t3_decomp/build_two_knockout.py (the captain's order of 2026-11-26: T3's
refinement round on C2's command, the T3 REF3 design as its cold challenge fixes it): on `--form c2-s2`, exactly the
five option strings of the builder's one literal table R76_TABLE, each only with `--captain-ruling R76`, with
`--den-ruling R76` for the three the row lifts the 5% check for and never for the two weight-only strings:

  T3_c2_s2_s085     `--form c2-s2 --size-scale 0.85`                        g = 1, s = 0.85   lifted
  T3_c2_s2_s05      `--form c2-s2 --size-scale 0.5`                         g = 1, s = 0.5    lifted
  T3_c2_s2_g15      `--form c2-s2 --gata4-exponent 1.5`                     g = 1.5, s = 1
  T3_c2_s2_g05      `--form c2-s2 --gata4-exponent 0.5`                     g = 0.5, s = 1
  T3_c2_s2_g2_s067  `--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67`   g = 2, s = 0.67   lifted

The table read off its strings; the factor (line 3's size-held product at g, value for value, then ** s) against its
formula and against C2's, line 3's and F2's arrays; the range guard and the 5% den check (exit 6, nothing written), the
latter gating the two weight-only files and lifted for the three; exit 5; the R76 check first and only with
`--captain-ruling R76`; every near-miss string, every other ruling and every den lift it may not take refused; the
ledger's row R76 holding the table (each string whole and paired with its file, one lifting part naming exactly the
three lifted strings and `--den-ruling R76`), on synthetic rows and on the repository's row R76 (its first cell read
by code, never printed); the files built on C2's cells, twice to one sha256; and every earlier command never touching
the line.

Offline, synthetic and fast: the throwaway world of tests/test_t3_two_knockout_build.py (a small synthetic WT file,
synthetic `transfer_mab21l2_all` and `gata4_pooled` tables, a decision file with the frozen hashes), runs whose two
tables are built here so that each check passes or fails by construction, and small ledger files under tmp_path.
Nothing here builds from released data or pins bytes built on one machine.
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
REPO_LEDGER = ROOT / "experiments" / "harness" / "TRIES.md"


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
SEED, N, G = th.SEED, th.N, th.G
PHI_M, PHI_G = th.PHI_M, th.PHI_G
INF = math.inf

TABLE = (
    ("T3_c2_s2_s085", "--form c2-s2 --size-scale 0.85", True),
    ("T3_c2_s2_s05", "--form c2-s2 --size-scale 0.5", True),
    ("T3_c2_s2_g15", "--form c2-s2 --gata4-exponent 1.5", False),
    ("T3_c2_s2_g05", "--form c2-s2 --gata4-exponent 0.5", False),
    ("T3_c2_s2_g2_s067", "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67", True),
)
NAMES = [n for n, _, _ in TABLE]
OPT = {n: o for n, o, _ in TABLE}
LIFTED = [n for n, _, lifted in TABLE if lifted]
WEIGHT_ONLY = [n for n, _, lifted in TABLE if not lifted]
KW = {"T3_c2_s2_s085": {"size_scale": 0.85}, "T3_c2_s2_s05": {"size_scale": 0.5},
      "T3_c2_s2_g15": {"gata4_exponent": 1.5}, "T3_c2_s2_g05": {"gata4_exponent": 0.5},
      "T3_c2_s2_g2_s067": {"gata4_exponent": 2.0, "size_scale": 0.67}}
POINTS = {"T3_c2_s2_s085": (1.0, 0.85), "T3_c2_s2_s05": (1.0, 0.5), "T3_c2_s2_g15": (1.5, 1.0),
          "T3_c2_s2_g05": (0.5, 1.0), "T3_c2_s2_g2_s067": (2.0, 0.67)}
LINE_FUNCS = ("r76_name", "r76_lift_part", "r76_row_holds_the_table", "r76_factor", "r76_guards", "run_r76_guards",
              "r76_loader", "build_r76")

ROW_R22 = ("| R22 (the captain's order of a day: a search; ruling, no run, not a try; the T3 two-knockout builder is "
           "to take `--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`; by T3's third search "
           "line also `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`; by its fourth line also "
           "`--form c2-s2 --gata4-exponent inf`, `--form c2-s2 --gata4-exponent inf --size-scale 0.67` and "
           "`--form c2-s2 --size-scale 0.67`, only with `--captain-ruling R22`) | T3 | x |\n")
ROW_R47 = ("| R47 (the captain's order of a day, ruling, no run: \"yes\"; for each of `--form c2-s2 --gata4-exponent "
           "inf`, `--form c2-s2 --gata4-exponent inf --size-scale 0.67` and `--form c2-s2 --size-scale 0.67`, the 5% "
           "check gates nothing for it; not a try) | T3 | y |\n")


def r76_row(names=tuple(NAMES), lifted=tuple(LIFTED), opener="R76 (the captain's order of a day",
            ruling="ruling, no run", captain="`--captain-ruling R76`", den="`--den-ruling R76`",
            words="the 5% check gates nothing for it", pairs=None) -> str:
    """A row R76 in the repository row's shape (the captain's words quoted first, then the strings with their files,
    then the one lifting part)."""
    pairs = pairs if pairs is not None else ", ".join(f"`{OPT[n]}` ({n})" for n in names)
    return (f"| {opener}, {ruling}: \"Yes to the five files under R76, sent at once: my named exception on T3's depth "
            "for each; for three my named exception for the Gata4 weight, as R41's; for three the 5% check gates "
            "nothing.\" (his pick in a dialog); T3's refinement round on C2's command; the T3 two-knockout builder is "
            f"to take exactly these strings, each only with {captain}: {pairs}; for each of "
            + ", ".join(f"`{OPT[n]}`" for n in lifted) + f", {words} ({den}); decision log a day (a); not a try) "
            "| T3 | z |\n")


ROW_R76 = r76_row()
ROWS = ROW_R22 + ROW_R47 + ROW_R76


def ledger(tmp_path: Path, rows: str = ROWS, name: str = "TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text(f"# a throwaway ledger\n| R21 (another row) | T2 | y |\n{rows}| R80 (a later row) | T3 | w |\n")
    return p


def repo_first_cell(row: str) -> str:
    """A row's first cell as the repository's ledger holds it, read by code (never printed)."""
    lines = [ln for ln in REPO_LEDGER.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    assert len(lines) == 1, f"the repository's ledger holds no single row {row}"
    return lines[0][2:].split("|", 1)[0]


def alternating(hi: float, lo: float, start: int, stop: int) -> np.ndarray:
    """1 on every gene but start..stop-1, where it alternates hi, lo, hi, lo."""
    v = np.ones(G)
    v[start:stop:2], v[start + 1:stop:2] = hi, lo
    return v


# MILD: every factor within about 2% of 1, so both checks pass for every file at any den. FAR: phi_G well above 1 on
# most genes, so each of the five points moves the factor's level away from C2's (each file's den-weighted total over
# 10% from C2's at the dens used here) while its range stays within C2's bounds (asserted where used); drawn once from a
# fixed seed, found by a search over seeds for that property.
MILD = (alternating(1.02, 1 / 1.02, 0, 10), alternating(1.02, 1 / 1.02, 10, 20))
_FAR_RNG = np.random.default_rng(284)
FAR = (np.exp(_FAR_RNG.normal(0.0, 0.5, G)), np.exp(_FAR_RNG.normal(0.8, 0.3, G)))
TABLES = {"mild": MILD, "far": FAR}


@pytest.fixture(scope="module")
def runs(world):
    """One run per table set (same genes, same WT record) and its two tables' sha256."""
    out = {}
    for key, (m, g) in TABLES.items():
        run = world.d / f"run_r76_{key}"
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


def build_r(world, tmp_path: Path, name: str, out: Path, rows: str = ROWS, **kw) -> dict:
    """Row R76's file `name` by its own command (with --den-ruling R76 where the row lifts it)."""
    lift = {"den_ruling": "R76"} if name in LIFTED else {}
    return build(world, out, captain_ruling="R76", tries=ledger(tmp_path, rows, name=f"{out.stem}_TRIES.md"),
                 **KW[name], **{**lift, **kw})


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


def ln_factor_by_hand(m, t, g: float, s: float = 1.0) -> np.ndarray:
    """The family's ln factor, written out: s c_g (ln phi_M + g ln phi_G), c_g = SD(ln M) / SD(ln M + g ln G)."""
    lm, lg = np.log(m), np.log(t)
    return s * (sd(lm) / sd(lm + g * lg)) * (lm + g * lg)


# ------------------------------------------------------------------------------------------------- the one table


def test_the_one_literal_table():
    assert bk.R76_TABLE == TABLE and (bk.R76_ROW, bk.R76_FORM, bk.R76_DEN_OPTION) == ("R76", "c2-s2",
                                                                                      "--den-ruling R76")
    assert bk.R76_FILES == {
        "T3_c2_s2_s085": ("--form c2-s2 --size-scale 0.85", None, 0.85, True),
        "T3_c2_s2_s05": ("--form c2-s2 --size-scale 0.5", None, 0.5, True),
        "T3_c2_s2_g15": ("--form c2-s2 --gata4-exponent 1.5", 1.5, None, False),
        "T3_c2_s2_g05": ("--form c2-s2 --gata4-exponent 0.5", 0.5, None, False),
        "T3_c2_s2_g2_s067": ("--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67", 2.0, 0.67, True)}
    assert list(bk.R76_FILES) == NAMES and len(set(OPT.values())) == 5
    # the design's struck s11 is not in it, nor any string an earlier row named on C2's command
    assert "T3_c2_s2_s11" not in bk.R76_FILES and not any("1.1" in o for o in OPT.values())
    earlier = [*bk.GONLY_OPTIONS.values(), bk.C2G_OPTION.strip("`"), bk.G2_OPTION.strip("`"),
               *bk.DEPTH_OPTIONS.values(), bk.DEPTH_FORM_OPTION, *(n.strip("`") for n in bk.R36_NAMES)]
    for a in OPT.values():
        for b in earlier:
            assert a != b and f"`{a}`" not in f"`{b}`" and f"`{b}`" not in f"`{a}`", (a, b)
    # nor C1 (g = 1 at s = 1/c) or C2 itself: every point moves g or s off C2's (1, 1)
    assert all(p != (1.0, 1.0) for p in POINTS.values())


@pytest.mark.parametrize("name,option,why", [
    ("T3_c2_s2_s085", "--form c2-s1 --size-scale 0.85", "not a row R76 option string"),
    ("T3_c2_s2_s085", "--size-scale 0.85", "not a row R76 option string"),
    ("T3_c2_s2_s085", "--form c2-s2", "not a row R76 option string"),
    ("T3_c2_s2_s085", "--form c2-s2 --size-scale", "not a row R76 option string"),
    ("T3_c2_s2_s085", "--form c2-s2 --size-scale .85", "not a row R76 option string"),
    ("T3_c2_s2_g2", "--form c2-s2 --gata4-exponent 2", "not a row R76 option string"),
    ("T3_c2_s2_g2_s067", "--form c2-s2 --size-scale 0.67 --gata4-exponent 2.0", "not a row R76 option string"),
    ("T3_c2_s2_g2", "--form c2-s2 --gata4-exponent 2.0 --gata4-exponent 2.0", "not a row R76 option string"),
    ("T3_c2_s2_x05", "--form c2-s2 --depth-exponent 0.5", "not a row R76 option string"),
    ("T3_c2_s2_s085", "--form c2-s2 --size-scale 0.85 --captain-ruling R76", "not a row R76 option string"),
    ("T3_c2_s2_s85", "--form c2-s2 --size-scale 0.85", "does not match its option string"),
    ("T3_c2_s2_g2_s067", "--form c2-s2 --gata4-exponent 2.0", "does not match its option string"),
    ("T3_c2_s2_g20_s067", "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67", "does not match"),
    ("T3_c2_s2_s05", "--form c2-s2 --size-scale 0.50", "does not match its option string"),
])
def test_the_table_is_read_off_its_strings_and_a_mismatch_stops_the_import(name, option, why):
    with pytest.raises(ValueError, match=why):
        bk._r76_values(name, option)


def test_the_tags_and_values():
    assert [bk._r76_tag(v) for v in ("2.0", "1.5", "0.5", "0.85", "0.67", "10.0")] == ["2", "15", "05", "085", "067",
                                                                                        "10"]
    for name, opt, _ in TABLE:
        g, s = bk._r76_values(name, opt)
        assert (1.0 if g is None else g, 1.0 if s is None else s) == POINTS[name]
        assert (g is None) == ("--gata4-exponent" not in opt) and (s is None) == ("--size-scale" not in opt)


# ------------------------------------------------------------------------------------------------------ the factor


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("name", NAMES)
def test_the_factor_is_line_3_s_size_held_product_then_the_size_scale(seed, name):
    rng = np.random.default_rng(seed)
    m = np.exp(rng.normal(0.05, 0.6, 300))
    t = np.exp(rng.normal(-0.02, 0.35, 300))
    g, s = POINTS[name]
    phi, rec = bk.r76_factor(m, t, g, s)
    held, hrec = bk.c2_gata4_weighted(m, t, g)
    assert phi.tobytes() == (held if s == 1.0 else held ** s).tobytes()  # the operations themselves, bit for bit
    assert np.allclose(np.log(phi), ln_factor_by_hand(m, t, g, s), rtol=1e-11, atol=1e-14)
    assert sd(np.log(phi)) == pytest.approx(s * sd(np.log(m)), rel=1e-12)  # held at C2's size, then times s
    assert rec["c"] == hrec["c"] == pytest.approx(sd(np.log(m)) / sd(np.log(m) + g * np.log(t)), rel=1e-12)
    assert rec["g"] == g and rec["size_scale"] == s
    assert rec["mab21l2_exponent"] == pytest.approx(s * rec["c"], rel=1e-15)
    assert rec["gata4_weight"] == pytest.approx(s * g * rec["c"], rel=1e-15)
    assert rec["min"] == float(phi.min()) and rec["max"] == float(phi.max())


def test_the_factor_at_c2_s_weight_is_c2_s_and_f2_s_and_at_2_line_3_s():
    """g = 1 gives combine('c2', 1.0)'s array (C2's) value for value; (1, 0.67) gives F2's array (the options on C2's
    command, c2_gonly_factor); (2, 1) gives line 3's (T3_c2_s2_g2's); the other points of c2_gonly_factor as they
    were."""
    rng = np.random.default_rng(7)
    for m, t in [(PHI_M, PHI_G), MILD, FAR, (np.exp(rng.normal(0, 0.6, 200)), np.exp(rng.normal(0, 0.4, 200)))]:
        c2, _ = bk.combine("c2", 1.0, m, t)
        assert bk.r76_factor(m, t, 1.0)[0].tobytes() == c2.tobytes()
        assert bk.r76_factor(m, t, 1.0, 0.67)[0].tobytes() == bk.c2_gonly_factor(m, t, 1.0, 0.67)[0].tobytes()
        assert bk.r76_factor(m, t, 1.0, 0.85)[0].tobytes() == (c2 ** 0.85).tobytes()
        assert bk.r76_factor(m, t, 2.0)[0].tobytes() == bk.c2_gata4_weighted(m, t, 2.0)[0].tobytes()
        x = bk.r76_factor(m, t, 2.0, 0.67)[0]
        assert x.tobytes() == (bk.c2_gata4_weighted(m, t, 2.0)[0] ** 0.67).tobytes()
    # the F-files' own function still takes only its points: R76's weights are not its
    for g in (1.5, 0.5, 2.0):
        with pytest.raises(bk.Refused, match="not this line's"):
            bk.c2_gonly_factor(PHI_M, PHI_G, g)


def test_the_factor_refuses_an_infinite_weight_bad_scales_and_other_genes():
    m = np.array([2.0, 0.5, 1.0])
    for g in (INF, -INF, float("nan")):
        with pytest.raises(bk.Refused, match="not finite"):
            bk.r76_factor(m, m[::-1], g)
    for s in (0.0, -0.5, INF, float("nan")):
        with pytest.raises(bk.Refused, match="size scale"):
            bk.r76_factor(m, m[::-1], 1.0, s)
    with pytest.raises(bk.Refused, match="same genes"):
        bk.r76_factor(m, np.array([2.0, 0.5]), 1.5)
    with pytest.raises(bk.Refused, match="finite and positive"):
        bk.r76_factor(m, np.array([2.0, 0.0, 1.0]), 1.5)


# ------------------------------------------------------------------------------------------------------ the checks


@pytest.mark.parametrize("key", list(TABLES))
@pytest.mark.parametrize("name", NAMES)
def test_the_two_checks_by_hand(key, name):
    m, t = TABLES[key]
    g, s = POINTS[name]
    for den in (np.ones(G), np.linspace(0.1, 3.0, G)):
        r = bk.r76_guards(m, t, g, s, den)
        phi = np.exp(ln_factor_by_hand(m, t, g, s))
        c2 = np.exp(ln_factor_by_hand(m, t, 1.0))
        rg, d = r["range"], r["den_check"]
        assert rg["factor_min"] == pytest.approx(phi.min(), rel=1e-11)
        assert rg["factor_max"] == pytest.approx(phi.max(), rel=1e-11)
        assert rg["bound_min"] == rg["c2_min"] / 2.0 and rg["bound_max"] == rg["c2_max"] * 2.0
        assert rg["passed"] is (rg["factor_min"] >= rg["bound_min"] and rg["factor_max"] <= rg["bound_max"])
        tot = den.sum()
        assert d["den_total_c2"] == pytest.approx(den @ c2 / tot, rel=1e-11)
        assert d["den_total_factor"] == pytest.approx(den @ phi / tot, rel=1e-11)
        assert d["den_shift_vs_c2"] == pytest.approx((den @ phi) / (den @ c2) - 1.0, rel=1e-10, abs=1e-14)
        assert d["within"] is (abs(d["den_shift_vs_c2"]) <= 0.05) and d["max_den_shift"] == 0.05
        assert d["passed"] is d["within"] and d["lifted_by"] is None and d["limit"]["status"] == "gates"
        assert r["passed"] is (rg["passed"] and d["within"]) and d["rule"] == bk.R76_DEN_RULE
        assert r["factor_sha256"] == hashlib.sha256(np.ascontiguousarray(bk.r76_factor(m, t, g, s)[0])
                                                    .tobytes()).hexdigest()
        assert set(r) == {"g", "size_scale", "c", "powers_applied", "range", "den_check", "passed", "factor_sha256"}
        lift = {"row": "R76", "option": f"`{OPT[name]}`"}
        lr = bk.r76_guards(m, t, g, s, den, lift)
        assert lr["den_check"]["passed"] is True and lr["den_check"]["within"] is d["within"]
        assert lr["den_check"]["den_shift_vs_c2"] == d["den_shift_vs_c2"] and lr["den_check"]["lifted_by"] == lift
        assert lr["den_check"]["limit"]["status"] == "lifted by row R76: the 5% check gates nothing for it"
        assert lr["passed"] is rg["passed"] and lr["range"] == rg  # the lift never touches the range guard
    if key == "mild":
        assert r["passed"]
    else:  # FAR: every point far from C2's level, its range within C2's bounds
        assert not r["den_check"]["within"] and r["range"]["passed"]


def test_the_range_bound_is_inclusive_and_a_lift_does_not_lift_it(monkeypatch):
    c2, _ = bk.combine("c2", 1.0, *MILD)
    for edge, ok in (("max_at", True), ("max_beyond", False), ("min_at", True), ("min_beyond", False)):
        f = np.ones(G)
        if edge.startswith("max"):
            f[0] = 2.0 * c2.max() if edge == "max_at" else np.nextafter(2.0 * c2.max(), INF)
        else:
            f[0] = c2.min() / 2.0 if edge == "min_at" else np.nextafter(c2.min() / 2.0, 0.0)
        monkeypatch.setattr(bk, "r76_factor", lambda *a, _f=f, **k: (_f, {"g": 2.0, "size_scale": 0.67, "c": 1.0,
                                                                          "mab21l2_exponent": 0.67,
                                                                          "gata4_weight": 1.34}))
        r = bk.r76_guards(*MILD, 2.0, 0.67, np.ones(G), lift={"row": "R76"})
        assert r["range"]["passed"] is ok and r["passed"] is ok and r["den_check"]["passed"]


def test_a_failed_check_says_which():
    r = bk.r76_guards(*FAR, 1.5, 1.0, np.ones(G))
    says = bk.gonly_failure(r)
    assert says.startswith("the 5% den check: the file's den-weighted total is ") and says.endswith(
        "against C2's factor's, beyond 5%, and no --den-ruling row lifts it for this file: not built")


# ------------------------------------------------------------------------------------------------- the ruling row


def test_without_the_ruling_nothing_of_the_line_is_read(monkeypatch):
    monkeypatch.setattr(bk, "_outside_form", lambda: pytest.fail("read on a command without --captain-ruling R76"))
    for name in LINE_FUNCS:
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached without R76"))
    for form in [*bk.FORMS, "", None]:
        for ruling in (None, "R22", "R36", "R47", "R53", "R41", "r76", "R76 ", " R76", "R760", "R7", "R76b"):
            for kw in ({}, *KW.values(), {"den_ruling": "R76"}, {"gata4_exponent": 2.0, "size_scale": 0.67,
                                                                  "den_ruling": "R76"}):
                assert bk.check_r76(form, captain_ruling=ruling, **kw) is None


@pytest.mark.parametrize("name", NAMES)
def test_each_string_is_taken_with_its_ruling_and_its_lift(tmp_path, name):
    p = ledger(tmp_path)
    lift = {"den_ruling": "R76"} if name in LIFTED else {}
    rec = bk.check_r76("c2-s2", captain_ruling="R76", tries=p, **KW[name], **lift)
    sha = t3d.sha256(p)
    want_lift = {"row": "R76", "option": f"`{OPT[name]}`", "words": "the 5% check gates nothing for it",
                 "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha} if name in LIFTED else None
    assert rec == {"row": "R76", "name": name, "file": f"{name}.h5ad", "option": f"`{OPT[name]}`",
                   "g": POINTS[name][0], "s": POINTS[name][1], "den_ruling": want_lift,
                   "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha}


def test_the_repository_s_row_r76_holds_the_table_and_opens_each_file(tmp_path):
    """Row R76's first cell as the repository's ledger holds it (read by code, never printed), in a throwaway ledger:
    it holds the table, its one lifting part names exactly the three lifted strings, and each file is taken."""
    first = repo_first_cell("R76")
    assert bk.r76_row_holds_the_table(first)
    assert sorted(bk.r76_lift_part(first)) == sorted([*(OPT[n] for n in LIFTED), "--den-ruling R76"])
    p = ledger(tmp_path, f"| {first}| T3 | x |\n")
    for name in NAMES:
        lift = {"den_ruling": "R76"} if name in LIFTED else {}
        assert bk.check_r76("c2-s2", captain_ruling="R76", tries=p, **KW[name], **lift)["name"] == name


def test_the_repository_s_rows_still_open_the_identity_rebuilds_commands(tmp_path):
    """The launcher's identity rebuilds on the repository's rows R22 and R47 (first cells read by code): line 3's
    T3_c2_s2_g2, and F2 and F1 with --den-ruling R47."""
    p = ledger(tmp_path, "".join(f"| {repo_first_cell(r)}| T3 | x |\n" for r in ("R22", "R47")))
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", p)["option"] == "`--form c2-s2 --gata4-exponent 2.0`"
    for kw, name in (({"size_scale": 0.67}, "F2"), ({"gata4_exponent": INF}, "F1")):
        rec = bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, den_ruling="R47", **kw)
        assert rec["name"] == name and rec["den_ruling"]["row"] == "R47"


@pytest.mark.parametrize("g,s", [
    (None, None), (1.0, 0.85), (1.5, 1.0), (0.5, 1.0), (1.0, None), (None, 1.0), (1.0, 1.0),
    (None, 0.67), (2.0, None), (INF, None), (INF, 0.67), (INF, 0.85), (None, 0.850001), (1.50001, None),
    (0.5, 0.5), (1.5, 0.67), (2.0, 0.5), (2.0, 0.85), (2.0, 0.670001), (2.000001, 0.67), (0.85, None), (None, 1.5),
    (-0.5, None), (None, -0.5), (None, 1.1), (3.0, None), (0.0, None), (float("nan"), None), (None, float("nan")),
    (True, None), (None, True), ("x", None), (None, "0,85"), (None, ""),
])
def test_a_near_miss_string_is_refused(tmp_path, g, s):
    """Every string near one of the five (a value off, an option added, dropped or swapped, line 3's and F2's own
    strings, the struck s11) is refused, whatever --den-ruling says."""
    p = ledger(tmp_path)
    for den in (None, "R76"):
        with pytest.raises(bk.Refused, match="(not one of row R76's five option strings|is refused under)"):
            bk.check_r76("c2-s2", g, s, "R76", p, den_ruling=den)


@pytest.mark.parametrize("form", [f for f in bk.FORMS if f != "c2-s2"])
def test_another_form_is_refused(tmp_path, form):
    for name in NAMES:
        with pytest.raises(bk.Refused, match="on --form c2-s2 only"):
            bk.check_r76(form, captain_ruling="R76", tries=ledger(tmp_path), **KW[name])


@pytest.mark.parametrize("extra", [{"composition_from": "allowed-s1"}, {"depth_from": "allowed-s1"},
                                   {"depth_exponent": 0.5}, {"depth_form": "global"}])
def test_the_composition_and_depth_options_beside_the_ruling_are_refused(tmp_path, extra):
    for name in NAMES:
        lift = {"den_ruling": "R76"} if name in LIFTED else {}
        with pytest.raises(bk.Refused, match="--captain-ruling R76 takes no --(composition-from|depth-)"):
            bk.check_r76("c2-s2", captain_ruling="R76", tries=ledger(tmp_path), **KW[name], **lift, **extra)


@pytest.mark.parametrize("den", [None, "R47", "R22", "R30", "r76", "R76 ", "R760", "", "R7"])
@pytest.mark.parametrize("name", LIFTED)
def test_a_lifted_string_is_taken_only_with_den_ruling_r76(tmp_path, name, den):
    with pytest.raises(bk.Refused, match="is taken only with --den-ruling R76"):
        bk.check_r76("c2-s2", captain_ruling="R76", tries=ledger(tmp_path), den_ruling=den, **KW[name])


@pytest.mark.parametrize("den", ["R76", "R47", "R22", "r76", ""])
@pytest.mark.parametrize("name", WEIGHT_ONLY)
def test_the_den_lift_is_refused_for_the_two_weight_only_strings(tmp_path, name, den):
    """check_den_lift's reading (the words anywhere in the cell) would lift them under row R76: the table does not."""
    p = ledger(tmp_path)
    first = bk.ledger_first_cell("R76", p)
    assert bk.GONLY_LIFT_WORDS in first and bk.names_whole(first, OPT[name])  # the hazard the table closes
    with pytest.raises(bk.Refused, match="lifts the 5% check for three strings only, not this one"):
        bk.check_r76("c2-s2", captain_ruling="R76", tries=p, den_ruling=den, **KW[name])


ROW_REFUSALS = [
    ("", "no row R76"),
    (ROW_R76 * 2, "two rows R76"),
    (ROW_R76.replace("| R76 (", "| R76b ("), "another row id"),
    (r76_row(opener="R76 (the captain's ruling of a day"), "not the captain's order"),
    (r76_row(opener="R76 (the captain's named exception of a day"), "not the captain's order"),
    (r76_row(ruling="a run"), "no ruling, no run"),
    (r76_row(captain="--captain-ruling R76"), "the ruling without backticks"),
    (r76_row(captain="`--captain-ruling R22`"), "another ruling"),
    (r76_row(den="--den-ruling R76"), "the lift's ruling without backticks"),
    (r76_row(den="`--den-ruling R47`"), "another lifting row"),
    (r76_row(words="the 5% check gates nothing for them"), "the lifting words changed"),
    (r76_row(words="the 5% check gates nothing"), "the lifting words cut"),
    (r76_row(lifted=LIFTED[:2]), "the lifting part misses one lifted string"),
    (r76_row(lifted=(*LIFTED, "T3_c2_s2_g15")), "the lifting part names a weight-only string"),
    (r76_row(lifted=NAMES), "the lifting part names every string"),
    (ROW_R76.replace(" (`--den-ruling R76`); decision", " (`--den-ruling R76`); the 5% check gates nothing for it; "
                     "decision"), "two parts say the words"),
    (r76_row(names=NAMES[:4]), "a string missing from the table"),
    (r76_row(names=NAMES[1:]), "the first string missing"),
    (r76_row(pairs=", ".join(f"`{OPT[n]}` ({m})" for n, m in zip(NAMES, [*NAMES[1:], NAMES[0]], strict=True))),
     "the strings paired with other files"),
    (r76_row(pairs=", ".join(f"`{OPT[n]}`" for n in NAMES)), "the strings without their files"),
    (r76_row(pairs=", ".join(f"{OPT[n]} ({n})" for n in NAMES)), "the strings without backticks"),
    (r76_row(pairs=", ".join(f"`a`{OPT[n]}`b` ({n})" for n in NAMES)), "not as whole strings"),
    (r76_row(pairs=", ".join(f"`{OPT[n]} ` ({n})" for n in NAMES)), "a trailing space inside the backticks"),
    (r76_row(pairs=", ".join(f"`{OPT[n].replace('c2-s2', 'c2-s1')}` ({n})" for n in NAMES)), "another form"),
    ("| R76 (the captain's order of a day, ruling, no run; not a try) | T3 | "
     + ", ".join(f"`{OPT[n]}` ({n})" for n in NAMES) + " `--captain-ruling R76`; for each of "
     + ", ".join(f"`{OPT[n]}`" for n in LIFTED) + ", the 5% check gates nothing for it (`--den-ruling R76`) |\n",
     "the table outside the first cell"),
]


@pytest.mark.parametrize("rows,why", ROW_REFUSALS, ids=[w for _, w in ROW_REFUSALS])
def test_a_row_r76_that_does_not_hold_the_table_opens_no_file(tmp_path, rows, why):
    p = ledger(tmp_path, ROW_R22 + ROW_R47 + rows)
    for name in NAMES:
        lift = {"den_ruling": "R76"} if name in LIFTED else {}
        with pytest.raises(bk.Refused, match="has no single row R76 whose first cell is the captain's order"):
            bk.check_r76("c2-s2", captain_ruling="R76", tries=p, **KW[name], **lift)


def test_the_row_check_is_build_outside_form_s_plus_the_table(tmp_path, monkeypatch):
    bo = bk._outside_form()
    seen = []
    real = bo.captain_order_names
    monkeypatch.setattr(bo, "captain_order_names", lambda row, option, path=None: seen.append((row, option, path))
                        or real(row, option, path))
    p = ledger(tmp_path)
    for name in NAMES:
        seen.clear()
        lift = {"den_ruling": "R76"} if name in LIFTED else {}
        assert bk.check_r76("c2-s2", captain_ruling="R76", tries=p, **KW[name], **lift)["name"] == name
        assert seen == [("R76", f"`{OPT[name]}`", p)]


@pytest.mark.parametrize("ruling", [None, "R22", "R47", "R41", "R53", "R36", "R30", "r76", "R76 "])
@pytest.mark.parametrize("name", NAMES)
def test_each_string_under_any_other_ruling_is_refused_as_before(world, tmp_path, monkeypatch, ruling, name):
    """Under --captain-ruling R22 only 0.67 is a size scale and only 2.0 (or inf) a weight on c2-s2; under any other
    ruling every line refuses; R76's own code is never reached."""
    for fn in LINE_FUNCS:
        monkeypatch.setattr(bk, fn, lambda *a, _n=fn, **k: pytest.fail(f"{_n} reached under {ruling!r}"))
    lift = {"den_ruling": "R76"} if name in LIFTED else {}
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused):
        build(world, out, captain_ruling=ruling, tries=ledger(tmp_path), **KW[name], **lift)
    assert nothing_written(out) and loader_restored()
    if ruling == "R22":
        says = "only 2.0" if "size_scale" not in KW[name] else ("only 0.67" if "gata4_exponent" not in KW[name]
                                                                 else "taken alone")
        with pytest.raises(bk.Refused, match=says):
            build(world, out, captain_ruling="R22", tries=ledger(tmp_path), **KW[name])


def test_line_3_s_and_f2_s_strings_under_r76_are_refused_and_build_as_before_under_r22(world, runs, tmp_path,
                                                                                      monkeypatch):
    for kw in ({"gata4_exponent": 2.0}, {"size_scale": 0.67}, {"gata4_exponent": INF}, {}):
        for den in (None, "R76"):
            out = tmp_path / "x.h5ad"
            with pytest.raises(bk.Refused, match="not one of row R76's five option strings"):
                build(world, out, captain_ruling="R76", tries=ledger(tmp_path), den_ruling=den, **kw)
            assert nothing_written(out)
    r = use(runs, monkeypatch, "mild")
    side = build(world, tmp_path / "f2.h5ad", run_dir=r.run, captain_ruling="R22", tries=ledger(tmp_path),
                 size_scale=0.67, den_ruling="R47")
    assert side["combined_build"]["c2_gonly_line"]["file"] == "T3_c2_s2_s067.h5ad"
    assert "r76_line" not in side["combined_build"]


# ------------------------------------------------------------------------------------------- the files at the build


@pytest.mark.parametrize("name", NAMES)
def test_each_file_is_c2_s_command_with_its_factor_and_records_both_checks(world, runs, tmp_path, monkeypatch, name):
    r = use(runs, monkeypatch, "mild")
    spy = Spy(monkeypatch)
    c2 = build(world, tmp_path / "c2.h5ad", run_dir=r.run)
    out = tmp_path / f"{name}.h5ad"
    side = build_r(world, tmp_path, name, out, run_dir=r.run)
    assert len(spy.seen) == 2
    want, rec = bk.r76_factor(r.m, r.g, *POINTS[name])
    assert np.array_equal(spy.seen[1], want) and np.array_equal(spy.seen[0], bk.combine("c2", 1.0, r.m, r.g)[0])
    # C2's command with its own option string: spec, inputs, streams, composition, seed, count and board are C2's
    assert side["spec"] == c2["spec"] and side["inputs"] == c2["inputs"] and side["spec_source"] == c2["spec_source"]
    assert side["inputs"]["rng_streams"] == ["submission:barred:s2", "submission:barred:s2:assemble"]
    assert side["construction"]["composition"] == c2["construction"]["composition"]
    assert side["construction"]["factor_exponent"] == 1.0 and side["output"]["check"]["ok"]
    assert side["construction"]["depth_ratios"] is None and set(side["spec"]) == {"factor", "composition"}
    cb, cb2 = side["combined_build"], c2["combined_build"]
    assert {k: v for k, v in cb.items() if k not in ("combination", "r76_line")} == {
        k: v for k, v in cb2.items() if k != "combination"}
    assert cb["form"] == "c2-s2" and cb["combination"]["c"] == rec["c"] and cb["depth_step"] is None
    line = cb["r76_line"]
    assert line["captain_ruling"]["row"] == "R76" and line["captain_ruling"]["name"] == name
    assert line["file"] == f"{name}.h5ad" and line["option"] == f"`{OPT[name]}`" and line["line"] == bk.R76_LINE
    assert (line["g"], line["size_scale"]) == POINTS[name] and line["c"] == rec["c"]
    assert line["powers_applied"] == {bk.TABLE_M: rec["mab21l2_exponent"], bk.TABLE_G: rec["gata4_weight"]}
    assert line["closed_lever"] == bk.R76_CLOSED_LEVER and line["list_entry"] == bk.R76_LIST_ENTRY
    assert line["c2_reference_sha256"] == "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"
    sha = hashlib.sha256(np.ascontiguousarray(want).tobytes()).hexdigest()
    assert line["factor_sha256_applied"] == line["guards"]["factor_sha256"] == sha
    gd = line["guards"]
    assert "genes" not in gd and gd["n_genes"] == G and gd["tables_sha256"] == r.pins
    assert gd["passed"] and gd["range"]["passed"] and gd["den_check"]["passed"] and gd["den_check"]["within"]
    if name in LIFTED:
        assert line["den_ruling"]["row"] == "R76" and gd["den_check"]["lifted_by"] == line["den_ruling"]
        assert gd["den_check"]["limit"]["status"] == "lifted by row R76: the 5% check gates nothing for it"
    else:
        assert line["den_ruling"] is None and gd["den_check"]["limit"]["status"] == "gates"
    want_g = bk.r76_guards(r.m, r.g, *POINTS[name], den_of(world), line["den_ruling"])

    def as_json(o):
        return json.loads(json.dumps(t3d.json_safe(o)))

    assert as_json({k: gd[k] for k in want_g}) == as_json(want_g)
    assert np.array_equal(np.array(gd["den"]), den_of(world)) and gd["wt_sha256"] == t3d.sha256(world.wt)
    on_disk = json.loads(out.with_suffix(".json").read_text())["combined_build"]["r76_line"]
    assert on_disk["file"] == f"{name}.h5ad" and on_disk["guards"]["range"] == gd["range"]
    assert json.loads(Path(f"{out}.stamp.json").read_text())["clean"] is True
    assert not any("cell" in k.lower() for k in line) and loader_restored()


@pytest.mark.parametrize("name", NAMES)
def test_each_file_draws_c2_s_cells_differs_in_expression_and_builds_one_sha256(world, runs, tmp_path, monkeypatch,
                                                                                 name):
    import anndata as ad

    r = use(runs, monkeypatch, "mild")
    build(world, tmp_path / "c2.h5ad", run_dir=r.run)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    sa = build_r(world, tmp_path, name, tmp_path / "a" / f"{name}.h5ad", run_dir=r.run)
    sb = build_r(world, tmp_path, name, tmp_path / "b" / f"{name}.h5ad", run_dir=r.run)
    assert (tmp_path / "a" / f"{name}.h5ad").read_bytes() == (tmp_path / "b" / f"{name}.h5ad").read_bytes()
    assert sa["output"]["sha256"] == sb["output"]["sha256"]
    a, c = ad.read_h5ad(tmp_path / "a" / f"{name}.h5ad"), ad.read_h5ad(tmp_path / "c2.h5ad")
    assert list(a.obs_names) == list(c.obs_names) and a.obs.equals(c.obs)
    assert np.array_equal(a.obsm["spatial_3D"], c.obsm["spatial_3D"])
    assert not np.array_equal(np.asarray(a.X), np.asarray(c.X))


def test_the_five_files_are_five_factors(world, runs, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    shas = {name: bk.factor_sha256(bk.r76_factor(r.m, r.g, *POINTS[name])[0]) for name in NAMES}
    c2 = bk.factor_sha256(bk.combine("c2", 1.0, r.m, r.g)[0])
    assert len(set(shas.values())) == 5 and c2 not in shas.values()


@pytest.mark.parametrize("name", NAMES)
def test_the_5_percent_check_stops_a_weight_only_file_and_gates_nothing_for_a_lifted_one(world, runs, tmp_path,
                                                                                         monkeypatch, name):
    """FAR: every file's den-weighted total is far from C2's and its range within C2's bounds. The two weight-only
    files stop (exit 6) before any cell is drawn, nothing written; the three lifted files build and record the
    figures, which gate nothing."""
    r = use(runs, monkeypatch, "far")
    pre = bk.r76_guards(r.m, r.g, *POINTS[name], den_of(world))
    assert not pre["den_check"]["within"] and pre["range"]["passed"]  # the premise
    out = tmp_path / f"{name}.h5ad"
    if name in WEIGHT_ONLY:
        monkeypatch.setattr(bk, "r76_loader", lambda *a, **k: pytest.fail("the build was reached past a failed check"))
        monkeypatch.setattr(bk.bs, "construct", lambda *a, **k: pytest.fail("a cell was drawn past a failed check"))
        with pytest.raises(bk.R76GuardFailed, match="the 5% den check"):
            build_r(world, tmp_path, name, out, run_dir=r.run)
        assert nothing_written(out) and loader_restored()
        assert issubclass(bk.R76GuardFailed, bk.GonlyGuardFailed)  # the command line's exit 6
    else:
        side = build_r(world, tmp_path, name, out, run_dir=r.run)
        d = side["combined_build"]["r76_line"]["guards"]["den_check"]
        assert d["passed"] and not d["within"] and abs(d["den_shift_vs_c2"]) > 0.05
        assert d["lifted_by"]["row"] == "R76" and side["output"]["check"]["ok"]


def test_the_range_guard_stops_a_lifted_file(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    c2, _ = bk.combine("c2", 1.0, r.m, r.g)
    wide = np.ones(G)
    wide[0] = 3.0 * c2.max()
    monkeypatch.setattr(bk, "r76_factor", lambda *a, **k: (wide, {"g": 2.0, "size_scale": 0.67, "c": 1.0,
                                                                  "mab21l2_exponent": 0.67, "gata4_weight": 1.34}))
    out = tmp_path / "T3_c2_s2_g2_s067.h5ad"
    with pytest.raises(bk.R76GuardFailed, match="the range guard"):
        build_r(world, tmp_path, "T3_c2_s2_g2_s067", out, run_dir=r.run)
    assert nothing_written(out) and loader_restored()


def test_a_factor_other_than_the_checks_exits_5_and_writes_nothing(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    real = bk.run_r76_guards
    monkeypatch.setattr(bk, "run_r76_guards", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.FactorMismatch, match="not the one the guards describe"):
        build_r(world, tmp_path, "T3_c2_s2_g15", out, run_dir=r.run)
    assert nothing_written(out)
    monkeypatch.setattr(bk, "run_r76_guards", lambda *a, **k: real(*a, **k) | {"genes": ["Z"]})
    with pytest.raises(bk.FactorMismatch, match="genes"):
        build_r(world, tmp_path, "T3_c2_s2_s05", out, run_dir=r.run)
    assert nothing_written(out) and loader_restored()


def test_the_checks_read_the_pinned_tables_only(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    monkeypatch.setattr(bk, "TABLE_SHA256", {**r.pins, bk.TABLE_G: "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match="pinned table"):
        build_r(world, tmp_path, "T3_c2_s2_g05", out, run_dir=r.run)
    assert nothing_written(out)


def test_an_existing_output_is_refused_before_the_checks(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    out = tmp_path / "T3_c2_s2_s085.h5ad"
    out.write_bytes(b"x")
    monkeypatch.setattr(bk, "run_r76_guards", lambda *a, **k: pytest.fail("the checks ran past an existing file"))
    with pytest.raises(bk.Refused, match="exists; pass --overwrite"):
        build_r(world, tmp_path, "T3_c2_s2_s085", out, run_dir=r.run)
    assert out.read_bytes() == b"x"


def test_an_unclean_tree_is_refused(world, tmp_path, own_repo):
    (own_repo / "stray.py").write_text("x = 1\n")
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match="not clean"):
        build_r(world, tmp_path, "T3_c2_s2_g15", out)
    assert nothing_written(out)


# ------------------------------------------------------------------------------- every earlier command as it was


OLD_SETS = [(f, {}, "world") for f in bk.FORMS] + [
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "world"),
    ("c2-s2", {"size_scale": 0.67, "captain_ruling": "R22", "den_ruling": "R47"}, "mild"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R22"}, "mild"),
    ("c2-s2", {"gata4_exponent": INF, "size_scale": 0.67, "captain_ruling": "R22"}, "mild"),
]


@pytest.mark.parametrize("form,kw,tables", OLD_SETS, ids=[f"{f}-{i}" for i, (f, _, _) in enumerate(OLD_SETS)])
def test_every_earlier_command_never_touches_the_line(world, runs, tmp_path, monkeypatch, form, kw, tables):
    run = world.run if tables == "world" else use(runs, monkeypatch, tables).run
    t = ledger(tmp_path)
    extra = dict(kw, tries=t) if kw else {}
    plain = build(world, tmp_path / "plain.h5ad", form, run_dir=run, **extra)
    for name in LINE_FUNCS:
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on {form}"))
    trapped = build(world, tmp_path / "trapped.h5ad", form, run_dir=run, **extra)
    assert (tmp_path / "plain.h5ad").read_bytes() == (tmp_path / "trapped.h5ad").read_bytes()
    assert plain["combined_build"] == trapped["combined_build"] and "r76_line" not in plain["combined_build"]


def test_line_3_s_file_never_touches_the_line(world, tmp_path, monkeypatch):
    """`--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22` on a run that passes line 3's own 5% guard."""
    run = world.d / "run_quarter_r76"
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
    assert plain["combined_build"] == trapped["combined_build"] and "r76_line" not in trapped["combined_build"]


# --------------------------------------------------------------------------------------------------------- CLI


def test_the_cli(world, runs, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    boards = th.cli_boards(tmp_path)
    r = use(runs, monkeypatch, "mild")

    def args(run, *extra, out):
        return ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
                str(boards), "--run-dir", str(run), "--seed", str(SEED), "--threads", "2", *extra, "--out", str(out)]

    opts = {n: OPT[n].split()[2:] for n in NAMES}
    runner = CliRunner()
    refusals = [[*OPT[n].split(), "--captain-ruling", "R22"] for n in NAMES]  # each string under R22: as before
    refusals += [[*OPT[n].split(), "--captain-ruling", "R76"] for n in LIFTED]  # no --den-ruling R76
    refusals += [[*OPT[n].split(), "--captain-ruling", "R76", "--den-ruling", "R76"] for n in WEIGHT_ONLY]
    refusals += [[*OPT[n].split(), "--captain-ruling", "R76", "--den-ruling", "R47"] for n in NAMES]
    refusals += [["--form", "c2-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R76"],
                 ["--form", "c2-s2", "--size-scale", "0.67", "--captain-ruling", "R76", "--den-ruling", "R76"],
                 ["--form", "c2-s2", "--captain-ruling", "R76"],
                 ["--form", "c2-s2", "--gata4-exponent", "1.0", "--size-scale", "0.85", "--captain-ruling", "R76",
                  "--den-ruling", "R76"],
                 ["--form", "c2-s2", "--size-scale", "1.1", "--captain-ruling", "R76", "--den-ruling", "R76"],
                 ["--form", "c1-s2", "--gata4-exponent", "1.5", "--captain-ruling", "R76"],
                 ["--form", "c2-s2", "--gata4-exponent", "1.5", "--captain-ruling", "R76", "--depth-from",
                  "allowed-s1"],
                 ["--form", "c2-s2", "--gata4-exponent", "0.5", "--captain-ruling", "R76", "--composition-from",
                  "allowed-s1"],
                 ["--form", "c2-s2", "--gata4-exponent", "1.5"], ["--form", "c2-s2", "--size-scale", "0.85"]]
    for extra in refusals:
        out = tmp_path / "r.h5ad"
        res = runner.invoke(bk.app, args(r.run, *extra, out=out))
        assert res.exit_code == 2 and "refused" in res.output, (extra, res.output)
        assert nothing_written(out)
    for name in NAMES:
        ok = tmp_path / f"{name}.h5ad"
        lift = ["--den-ruling", "R76"] if name in LIFTED else []
        res = runner.invoke(bk.app, args(r.run, "--form", "c2-s2", *opts[name], "--captain-ruling", "R76", *lift,
                                         out=ok))
        assert res.exit_code == 0, res.output
        den = "gates nothing for it (row R76)" if name in LIFTED else "passed"
        assert f"row R76's refinement round on C2's command: {name}.h5ad" in res.output
        assert f"the range guard passed; the 5% den check {den}" in res.output and res.output.splitlines()[-1] == (
            f"OK: {ok}")
        line = json.loads(ok.with_suffix(".json").read_text())["combined_build"]["r76_line"]
        assert line["file"] == f"{name}.h5ad" and line["captain_ruling"]["name"] == name
    # a check fails: exit 6, nothing written (FAR: the 5% den check of a weight-only file)
    use(runs, monkeypatch, "far")
    bad = tmp_path / "stop.h5ad"
    res = runner.invoke(bk.app, args(runs["far"].run, "--form", "c2-s2", "--gata4-exponent", "1.5",
                                     "--captain-ruling", "R76", out=bad))
    assert res.exit_code == 6 and "refused (the range guard or the 5% den check" in res.output, res.output
    assert "the 5% den check" in res.output and nothing_written(bad)
    # the factor is not the checks': exit 5, nothing written
    use(runs, monkeypatch, "mild")
    real = bk.run_r76_guards
    monkeypatch.setattr(bk, "run_r76_guards", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    mis = tmp_path / "mismatch.h5ad"
    res = runner.invoke(bk.app, args(r.run, "--form", "c2-s2", "--gata4-exponent", "0.5", "--captain-ruling", "R76",
                                     out=mis))
    assert res.exit_code == 5 and nothing_written(mis)
    monkeypatch.setattr(bk, "run_r76_guards", real)
    # a ledger without row R76: refused
    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path, ROW_R22 + ROW_R47, name="no_r76.md"))
    res = runner.invoke(bk.app, args(r.run, "--form", "c2-s2", "--gata4-exponent", "1.5", "--captain-ruling", "R76",
                                     out=tmp_path / "n.h5ad"))
    assert res.exit_code == 2 and "no single row R76" in res.output and nothing_written(tmp_path / "n.h5ad")
