"""Row R79's five files in experiments/t3_decomp/build_two_knockout.py (the captain's order of 2026-11-26: T3's
refinement round 4 on C2's command, the T3 REF4 design with its cold challenge's fixes, all five sent at once): on
`--form c2-s2`, exactly the five option strings of the builder's one literal table R79_TABLE, each only with
`--captain-ruling R79` and `--den-ruling R79` (the row lifts the 5% check for every one: each carries a size scale):

  T3_c2_s2_s035       `--form c2-s2 --size-scale 0.35`                       g = 1, s = 0.35
  T3_c2_s2_g05_s05    `--form c2-s2 --gata4-exponent 0.5 --size-scale 0.5`   g = 0.5, s = 0.5
  T3_c2_s2_gonly_s05  `--form c2-s2 --gata4-exponent inf --size-scale 0.5`   g = inf, s = 0.5
  T3_c2_s2_g2_s05     `--form c2-s2 --gata4-exponent 2.0 --size-scale 0.5`   g = 2, s = 0.5
  T3_c2_s2_g15_s05    `--form c2-s2 --gata4-exponent 1.5 --size-scale 0.5`   g = 1.5, s = 0.5

The table read off its strings by row R79's parser (inf exactly, tagged gonly; near misses such as 0.350 or Inf stop the
import); the factor at a finite g row R76's r76_factor value for value, at g = inf the options' c2_gonly_factor value
for value, each against its formula; the range guard (exit 6, nothing written) and the 5% den check, recorded and lifted
for all five; exit 5; the R79 check first and only with `--captain-ruling R79`; every near-miss value, every other
ruling, every other den ruling, the composition and depth options, and row R76's strings under R79 refused; the ledger's
row R79 holding the table (each string whole and paired with its file, one lifting part naming exactly the five strings
and `--den-ruling R79`), on synthetic rows and on the repository's row R79 (its first cell read by code, never printed),
with the lift-phrase trap: the row's quoted words say "the 5% check gates nothing." and must not count as a lifting
part; the files built on C2's cells, twice to one sha256, through the typer command line with the summary printing g =
inf; and every earlier command never touching the line.

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
import re
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
    ("T3_c2_s2_s035", "--form c2-s2 --size-scale 0.35", True),
    ("T3_c2_s2_g05_s05", "--form c2-s2 --gata4-exponent 0.5 --size-scale 0.5", True),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --gata4-exponent inf --size-scale 0.5", True),
    ("T3_c2_s2_g2_s05", "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.5", True),
    ("T3_c2_s2_g15_s05", "--form c2-s2 --gata4-exponent 1.5 --size-scale 0.5", True),
)
NAMES = [n for n, _, _ in TABLE]
OPT = {n: o for n, o, _ in TABLE}
KW = {"T3_c2_s2_s035": {"size_scale": 0.35}, "T3_c2_s2_g05_s05": {"gata4_exponent": 0.5, "size_scale": 0.5},
      "T3_c2_s2_gonly_s05": {"gata4_exponent": INF, "size_scale": 0.5},
      "T3_c2_s2_g2_s05": {"gata4_exponent": 2.0, "size_scale": 0.5},
      "T3_c2_s2_g15_s05": {"gata4_exponent": 1.5, "size_scale": 0.5}}
POINTS = {"T3_c2_s2_s035": (1.0, 0.35), "T3_c2_s2_g05_s05": (0.5, 0.5), "T3_c2_s2_gonly_s05": (INF, 0.5),
          "T3_c2_s2_g2_s05": (2.0, 0.5), "T3_c2_s2_g15_s05": (1.5, 0.5)}
FINITE = [n for n in NAMES if POINTS[n][0] != INF]
LINE_FUNCS = ("r79_name", "r79_row_holds_the_table", "r79_point", "r79_functions", "r79_factor", "build_r79")
R76_OPT = {"T3_c2_s2_s085": "--form c2-s2 --size-scale 0.85", "T3_c2_s2_s05": "--form c2-s2 --size-scale 0.5",
           "T3_c2_s2_g15": "--form c2-s2 --gata4-exponent 1.5", "T3_c2_s2_g05": "--form c2-s2 --gata4-exponent 0.5",
           "T3_c2_s2_g2_s067": "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67"}
R76_KW = {"T3_c2_s2_s085": {"size_scale": 0.85}, "T3_c2_s2_s05": {"size_scale": 0.5},
          "T3_c2_s2_g15": {"gata4_exponent": 1.5}, "T3_c2_s2_g05": {"gata4_exponent": 0.5},
          "T3_c2_s2_g2_s067": {"gata4_exponent": 2.0, "size_scale": 0.67}}
R76_LIFTED = ("T3_c2_s2_s085", "T3_c2_s2_s05", "T3_c2_s2_g2_s067")
QUOTE = ("\"Yes to the five files under R79, sent at once: my named exception on T3's depth for each; for three my "
         "named exception for the Gata4 weight, as R41's for T3_c2_s2_g2; for one my named exception for the Gata4 "
         "weight, as R47's for T3_c2_s2_gonly; the 5% check gates nothing.\"")

ROW_R22 = ("| R22 (the captain's order of a day: a search; ruling, no run, not a try; the T3 two-knockout builder is "
           "to take `--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`; by T3's third search "
           "line also `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`; by its fourth line also "
           "`--form c2-s2 --gata4-exponent inf`, `--form c2-s2 --gata4-exponent inf --size-scale 0.67` and "
           "`--form c2-s2 --size-scale 0.67`, only with `--captain-ruling R22`) | T3 | x |\n")
ROW_R47 = ("| R47 (the captain's order of a day, ruling, no run: \"yes\"; for each of `--form c2-s2 --gata4-exponent "
           "inf`, `--form c2-s2 --gata4-exponent inf --size-scale 0.67` and `--form c2-s2 --size-scale 0.67`, the 5% "
           "check gates nothing for it; not a try) | T3 | y |\n")
ROW_R76 = ("| R76 (the captain's order of a day, ruling, no run: \"Yes to five files under R76\"; T3's refinement "
           "round on C2's command; the T3 two-knockout builder is to take exactly these strings, each only with "
           "`--captain-ruling R76`: " + ", ".join(f"`{o}` ({n})" for n, o in R76_OPT.items()) + "; for each of "
           + ", ".join(f"`{R76_OPT[n]}`" for n in R76_LIFTED) + ", the 5% check gates nothing for it "
           "(`--den-ruling R76`); not a try) | T3 | w |\n")


def r79_row(names=tuple(NAMES), lifted=tuple(NAMES), opener="R79 (the captain's order of a day",
            ruling="ruling, no run", captain="`--captain-ruling R79`", den="`--den-ruling R79`",
            words="the 5% check gates nothing for it", pairs=None, quote=QUOTE) -> str:
    """A row R79 in the repository row's shape (the captain's words quoted first, ending "the 5% check gates
    nothing.", then the strings with their files, then the one lifting part)."""
    pairs = pairs if pairs is not None else ", ".join(f"`{OPT[n]}` ({n})" for n in names)
    return (f"| {opener}, {ruling}: {quote} (their pick in a dialog; the option recorded the words above as the "
            "captain's); T3's refinement round 4 on C2's command; the T3 two-knockout builder is to take exactly these "
            f"strings, each only with {captain}: {pairs}; for each of " + ", ".join(f"`{OPT[n]}`" for n in lifted)
            + f", {words} ({den}); decision log a day (h); not a try) | T3 | z |\n")


ROW_R79 = r79_row()
ROWS = ROW_R22 + ROW_R47 + ROW_R76 + ROW_R79


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
# 5% from C2's at the dens used here) while its range stays within C2's bounds (asserted where used); row R76's tests'
# draw, from the same fixed seed.
MILD = (alternating(1.02, 1 / 1.02, 0, 10), alternating(1.02, 1 / 1.02, 10, 20))
_FAR_RNG = np.random.default_rng(284)
FAR = (np.exp(_FAR_RNG.normal(0.0, 0.5, G)), np.exp(_FAR_RNG.normal(0.8, 0.3, G)))
TABLES = {"mild": MILD, "far": FAR}


@pytest.fixture(scope="module")
def runs(world):
    """One run per table set (same genes, same WT record) and its two tables' sha256."""
    out = {}
    for key, (m, g) in TABLES.items():
        run = world.d / f"run_r79_{key}"
        if not run.exists():
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
    """Row R79's file `name` by its own command (with --den-ruling R79, which the row gives every string)."""
    return build(world, out, captain_ruling="R79", tries=ledger(tmp_path, rows, name=f"{out.stem}_TRIES.md"),
                 **KW[name], **{"den_ruling": "R79", **kw})


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
    """The family's ln factor, written out: s c_g (ln phi_M + g ln phi_G), c_g = SD(ln M) / SD(ln M + g ln G); at
    g = inf, s (SD(ln M) / SD(ln G)) ln G (phi_M dropped)."""
    lm, lg = np.log(m), np.log(t)
    if g == INF:
        return s * (sd(lm) / sd(lg)) * lg
    return s * (sd(lm) / sd(lm + g * lg)) * (lm + g * lg)


def as_json(o):
    return json.loads(json.dumps(t3d.json_safe(o)))


# ------------------------------------------------------------------------------------------------- the one table


def test_the_one_literal_table_is_the_row_s():
    assert bk.R79_TABLE == TABLE and (bk.R79_ROW, bk.R79_FORM, bk.R79_DEN_OPTION) == ("R79", "c2-s2",
                                                                                      "--den-ruling R79")
    assert bk.R79_FILES == {
        "T3_c2_s2_s035": ("--form c2-s2 --size-scale 0.35", None, 0.35, True),
        "T3_c2_s2_g05_s05": ("--form c2-s2 --gata4-exponent 0.5 --size-scale 0.5", 0.5, 0.5, True),
        "T3_c2_s2_gonly_s05": ("--form c2-s2 --gata4-exponent inf --size-scale 0.5", INF, 0.5, True),
        "T3_c2_s2_g2_s05": ("--form c2-s2 --gata4-exponent 2.0 --size-scale 0.5", 2.0, 0.5, True),
        "T3_c2_s2_g15_s05": ("--form c2-s2 --gata4-exponent 1.5 --size-scale 0.5", 1.5, 0.5, True)}
    assert list(bk.R79_FILES) == NAMES and len(set(OPT.values())) == 5 and all(lift for *_, lift in bk.R79_TABLE)
    # no string any earlier row named, R76's five among them, contains or sits inside one of these
    earlier = [*bk.GONLY_OPTIONS.values(), bk.C2G_OPTION.strip("`"), bk.G2_OPTION.strip("`"),
               *bk.DEPTH_OPTIONS.values(), bk.DEPTH_FORM_OPTION, *(n.strip("`") for n in bk.R36_NAMES),
               *(o for _, o, _ in bk.R76_TABLE)]
    for a in OPT.values():
        for b in earlier:
            assert a != b and f"`{a}`" not in f"`{b}`" and f"`{b}`" not in f"`{a}`", (a, b)
    # every point moves s off 1, and none is C2's, C1's, an R76 point or an F-file's point
    r76_points = {(1.0 if g is None else g, 1.0 if s is None else s) for _, g, s, _ in bk.R76_FILES.values()}
    assert all(s != 1.0 for _, s in POINTS.values()) and not set(POINTS.values()) & (r76_points | {(1.0, 1.0)})
    assert not set(POINTS.values()) & set(bk.GONLY_POINTS.values())
    assert set(bk.R79_FILES).isdisjoint(bk.R76_FILES) and set(bk.R79_FILES).isdisjoint(
        n.removesuffix(".h5ad") for n in bk.GONLY_FILES.values())


@pytest.mark.parametrize("name,option,why", [
    # near misses of the strings: a spelling the row does not give stops the import
    ("T3_c2_s2_s035", "--form c2-s2 --size-scale 0.350", "does not match its option string"),
    ("T3_c2_s2_s035", "--form c2-s2 --size-scale .35", "not a row R79 option string"),
    ("T3_c2_s2_s035", "--form c2-s2 --size-scale 0.35 ", "not a row R79 option string"),
    ("T3_c2_s2_s035", "--form c2-s2  --size-scale 0.35", "not a row R79 option string"),
    ("T3_c2_s2_s035", " --form c2-s2 --size-scale 0.35", "not a row R79 option string"),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --gata4-exponent Inf --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --gata4-exponent INF --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --gata4-exponent infinity --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --gata4-exponent +inf --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --gata4-exponent -inf --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_gnan_s05", "--form c2-s2 --gata4-exponent nan --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_sonly", "--form c2-s2 --size-scale inf", "not a row R79 option string"),  # inf is a weight only
    ("T3_c2_s2_ginf_s05", "--form c2-s2 --gata4-exponent inf --size-scale 0.5", "does not match its option string"),
    ("T3_c2_s2_gonly", "--form c2-s2 --gata4-exponent inf --size-scale 0.5", "does not match its option string"),
    ("T3_c2_s2_gonly_s05", "--form c2-s2 --size-scale 0.5 --gata4-exponent inf", "not a row R79 option string"),
    ("T3_c2_s2_g2_s05", "--form c2-s2 --gata4-exponent 2 --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_g20_s05", "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.5", "does not match its option string"),
    ("T3_c2_s2_g15_s05", "--form c2-s2 --gata4-exponent 1.5 --size-scale 0.50", "does not match its option string"),
    ("T3_c2_s2_g15_s05", "--form c2-s1 --gata4-exponent 1.5 --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_g15_s05", "--gata4-exponent 1.5 --size-scale 0.5", "not a row R79 option string"),
    ("T3_c2_s2_g15_s05", "--form c2-s2 --gata4-exponent 1.5 --size-scale 0.5 --den-ruling R79", "not a row R79"),
    ("T3_c2_s2_s035", "--form c2-s2 --size-scale 0.35 --size-scale 0.35", "not a row R79 option string"),
    ("T3_c2_s2", "--form c2-s2", "not a row R79 option string"),
])
def test_the_table_is_read_off_its_strings_and_a_near_miss_stops_the_import(name, option, why):
    with pytest.raises(ValueError, match=why):
        bk._r79_values(name, option)


def test_the_tags_and_values():
    for name, opt, _ in TABLE:
        g, s = bk._r79_values(name, opt)
        assert (1.0 if g is None else g, s) == POINTS[name]
        assert (g is None) == ("--gata4-exponent" not in opt) and s is not None
    assert bk._r79_values("T3_c2_s2_gonly_s05", OPT["T3_c2_s2_gonly_s05"])[0] == INF
    assert (bk.R79_INF, bk.R79_INF_TAG) == ("inf", "only")
    # row R76's own parser is unchanged: it never took inf (the tag rule would have given ginf)
    with pytest.raises(ValueError, match="not a row R76 option string"):
        bk._r76_values("T3_c2_s2_gonly_s05", OPT["T3_c2_s2_gonly_s05"])
    assert bk._r76_tag("inf") == "inf"


# ------------------------------------------------------------------------------------------------------ the factor


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("name", NAMES)
def test_each_point_is_its_formula_and_its_branch_s_function_bit_for_bit(seed, name):
    rng = np.random.default_rng(seed)
    m = np.exp(rng.normal(0.05, 0.6, 300))
    t = np.exp(rng.normal(-0.02, 0.35, 300))
    g, s = POINTS[name]
    phi, rec = bk.r79_factor(m, t, g, s)
    assert np.allclose(phi, np.exp(ln_factor_by_hand(m, t, g, s)), rtol=1e-12, atol=0)
    assert sd(np.log(phi)) == pytest.approx(s * sd(np.log(m)), rel=1e-12)  # held at C2's size, then times s
    if g == INF:
        want, wrec = bk.c2_gonly_factor(m, t, INF, s)
        assert rec["g"] == "inf" and rec["mab21l2_exponent"] == 0.0
        assert rec["c"] == pytest.approx(sd(np.log(m)) / sd(np.log(t)), rel=1e-12)
    else:
        want, wrec = bk.r76_factor(m, t, g, s)
        assert rec["g"] == g and rec["mab21l2_exponent"] == pytest.approx(s * rec["c"], rel=1e-15)
        assert rec["c"] == pytest.approx(sd(np.log(m)) / sd(np.log(m) + g * np.log(t)), rel=1e-12)
    assert phi.tobytes() == want.tobytes() and rec == wrec  # the very function, value for value
    assert rec["size_scale"] == s and rec["min"] == float(phi.min()) and rec["max"] == float(phi.max())


def test_the_factor_at_r76_s_and_r47_s_points_is_their_arrays():
    """r79_factor at (1, 0.5) is R76's T3_c2_s2_s05 array, at (0.5, 1) R76's T3_c2_s2_g05 array, at (2, 0.67) R76's
    T3_c2_s2_g2_s067, at (inf, 0.67) R47's F3 (c2_gonly_factor) and at (inf, 1) F1; at (1, 1) C2's."""
    rng = np.random.default_rng(11)
    for m, t in [(PHI_M, PHI_G), MILD, FAR, (np.exp(rng.normal(0, 0.6, 200)), np.exp(rng.normal(0, 0.4, 200)))]:
        c2, _ = bk.combine("c2", 1.0, m, t)
        assert bk.r79_factor(m, t, 1.0, 1.0)[0].tobytes() == c2.tobytes()
        assert bk.r79_factor(m, t, 1.0, 0.5)[0].tobytes() == bk.r76_factor(m, t, 1.0, 0.5)[0].tobytes()
        assert bk.r79_factor(m, t, 0.5, 1.0)[0].tobytes() == bk.r76_factor(m, t, 0.5)[0].tobytes()
        assert bk.r79_factor(m, t, 2.0, 0.67)[0].tobytes() == bk.r76_factor(m, t, 2.0, 0.67)[0].tobytes()
        assert bk.r79_factor(m, t, INF, 0.67)[0].tobytes() == bk.c2_gonly_factor(m, t, INF, 0.67)[0].tobytes()
        assert bk.r79_factor(m, t, INF, 1.0)[0].tobytes() == bk.c2_gonly_factor(m, t, INF)[0].tobytes()
        assert bk.r79_factor(m, t, 1.0, 0.35)[0].tobytes() == (c2 ** 0.35).tobytes()
        assert bk.r79_factor(m, t, INF, 0.5)[0].tobytes() == (bk.c2_gonly_factor(m, t, INF)[0] ** 0.5).tobytes()


def test_the_branches_and_their_functions(monkeypatch):
    fin, inf = bk.r79_functions(1.5), bk.r79_functions(INF)
    assert (fin["factor"], fin["guards"], fin["run_guards"], fin["loader"]) == (
        bk.r76_factor, bk.r76_guards, bk.run_r76_guards, bk.r76_loader)
    assert (inf["factor"], inf["guards"], inf["run_guards"], inf["loader"]) == (
        bk.c2_gonly_factor, bk.gonly_guards, bk.run_gonly_guards, bk.gonly_loader)
    for g in (0.5, 1.0, 2.0, 1.5):
        assert bk.r79_functions(g)["branch"].startswith("finite g")
    assert bk.r79_functions(INF)["branch"].startswith("g = inf")
    # each is looked up when called: a replaced function is the one used
    monkeypatch.setattr(bk, "c2_gonly_factor", lambda *a, **k: ("x", {}))
    assert bk.r79_factor(PHI_M, PHI_G, INF, 0.5)[0] == "x"
    # R76's own function still refuses an infinite weight (the reason row R79 has its g = inf branch)
    with pytest.raises(bk.Refused, match="not finite"):
        bk.r76_factor(PHI_M, PHI_G, INF, 0.5)


# ------------------------------------------------------------------------------------------------------ the checks


@pytest.mark.parametrize("key", list(TABLES))
@pytest.mark.parametrize("name", NAMES)
def test_the_two_checks_by_hand_lifted_for_each(key, name):
    m, t = TABLES[key]
    g, s = POINTS[name]
    lift = {"row": "R79", "option": f"`{OPT[name]}`"}
    for den in (np.ones(G), np.linspace(0.1, 3.0, G)):
        r = bk.r79_functions(g)["guards"](m, t, g, s, den, lift)
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
        assert d["within"] is (abs(d["den_shift_vs_c2"]) <= 0.05) and d["passed"] is True and d["lifted_by"] == lift
        assert d["limit"]["status"] == "lifted by row R79: the 5% check gates nothing for it"
        assert r["passed"] is rg["passed"] and r["g"] == ("inf" if g == INF else g) and r["size_scale"] == s
        assert r["factor_sha256"] == hashlib.sha256(np.ascontiguousarray(bk.r79_factor(m, t, g, s)[0])
                                                    .tobytes()).hexdigest()
    if key == "mild":
        assert r["passed"] and d["within"]
    else:  # FAR: every point far from C2's level, its range within C2's bounds: lifted, so it builds
        assert not d["within"] and rg["passed"] and r["passed"]


# ------------------------------------------------------------------------------------------------- the ruling row


def test_without_the_ruling_nothing_of_the_line_is_read(monkeypatch):
    monkeypatch.setattr(bk, "_outside_form", lambda: pytest.fail("read on a command without --captain-ruling R79"))
    for name in LINE_FUNCS:
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached without R79"))
    for form in [*bk.FORMS, "", None]:
        for ruling in (None, "R22", "R36", "R47", "R53", "R41", "R76", "r79", "R79 ", " R79", "R790", "R7", "R79b"):
            for kw in ({}, *KW.values(), {"den_ruling": "R79"}, {"gata4_exponent": INF, "size_scale": 0.5,
                                                                  "den_ruling": "R79"}):
                assert bk.check_r79(form, captain_ruling=ruling, **kw) is None


@pytest.mark.parametrize("name", NAMES)
def test_each_string_is_taken_with_its_ruling_and_its_lift(tmp_path, name):
    p = ledger(tmp_path)
    rec = bk.check_r79("c2-s2", captain_ruling="R79", tries=p, den_ruling="R79", **KW[name])
    sha = t3d.sha256(p)
    g, s = POINTS[name]
    assert rec == {"row": "R79", "name": name, "file": f"{name}.h5ad", "option": f"`{OPT[name]}`",
                   "g": "inf" if g == INF else g, "s": s,
                   "den_ruling": {"row": "R79", "option": f"`{OPT[name]}`",
                                  "words": "the 5% check gates nothing for it",
                                  "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha},
                   "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha}
    assert bk.r79_point(rec) == POINTS[name]
    # the command line's inf, as typer reads it, and as a string the option would read
    if g == INF:
        for v in ("inf", float("inf"), np.float64("inf")):
            assert bk.check_r79("c2-s2", v, 0.5, "R79", p, den_ruling="R79")["name"] == name


def test_the_repository_s_row_r79_holds_the_table_and_opens_each_file(tmp_path):
    """Row R79's first cell as the repository's ledger holds it (read by code, never printed), in a throwaway ledger:
    it holds the table, its one lifting part names exactly the five strings and `--den-ruling R79`, and each file is
    taken."""
    first = repo_first_cell("R79")
    assert bk.r79_row_holds_the_table(first)
    assert sorted(bk.r76_lift_part(first)) == sorted([*OPT.values(), "--den-ruling R79"])
    # the table is the row's: its pairs of a whole string and a file, in the row's order, and no other
    assert re.findall(r"`([^`]*)` \((T3_[A-Za-z0-9_]+)\)", first) == [(o, n) for n, o, _ in bk.R79_TABLE]
    assert re.findall(r"`([^`]*)`", first).count("--captain-ruling R79") == 1
    p = ledger(tmp_path, f"| {first}| T3 | x |\n")
    for name in NAMES:
        assert bk.check_r79("c2-s2", captain_ruling="R79", tries=p, den_ruling="R79", **KW[name])["name"] == name


def test_the_lift_phrase_trap_on_the_repository_s_row_r79(tmp_path):
    """The trap the cold challenge found: the row quotes the captain's words, and a lifting part is a part between
    semicolons that says "the 5% check gates nothing for it". The real row's quoted words say "the 5% check gates
    nothing." (no "for it"): they are not a lifting part, and exactly one part lifts. Had the words said "... gates
    nothing for it" (the design's first wording), two parts would say it and every file would be refused."""
    first = repo_first_cell("R79")
    parts = first.split(";")
    loose = [p for p in parts if "the 5% check gates nothing" in p]
    strict = [p for p in parts if bk.GONLY_LIFT_WORDS in p]
    assert len(loose) == 2 and len(strict) == 1  # the quoted words' part is near the words, and does not count
    quoted = next(p for p in loose if p not in strict)
    assert "the 5% check gates nothing." in quoted and "`" not in quoted
    trapped = first.replace("the 5% check gates nothing.\"", "the 5% check gates nothing for it.\"", 1)
    assert trapped != first and bk.r76_lift_part(trapped) is None and not bk.r79_row_holds_the_table(trapped)
    p = ledger(tmp_path, f"| {trapped}| T3 | x |\n")
    for name in NAMES:
        with pytest.raises(bk.Refused, match="has no single row R79 whose first cell is the captain's order"):
            bk.check_r79("c2-s2", captain_ruling="R79", tries=p, den_ruling="R79", **KW[name])


def test_the_repository_s_rows_still_open_the_identity_rebuilds_commands(tmp_path):
    """The launcher's identity rebuilds on the repository's rows R22, R47 and R76 (first cells read by code): line 3's
    T3_c2_s2_g2, F2, F1 and F3 with --den-ruling R47, and R76's T3_c2_s2_s05 with --den-ruling R76."""
    p = ledger(tmp_path, "".join(f"| {repo_first_cell(r)}| T3 | x |\n" for r in ("R22", "R47", "R76")))
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", p)["option"] == "`--form c2-s2 --gata4-exponent 2.0`"
    for kw, name in (({"size_scale": 0.67}, "F2"), ({"gata4_exponent": INF}, "F1"),
                     ({"gata4_exponent": INF, "size_scale": 0.67}, "F3")):
        rec = bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, den_ruling="R47", **kw)
        assert rec["name"] == name and rec["den_ruling"]["row"] == "R47"
    rec = bk.check_r76("c2-s2", size_scale=0.5, captain_ruling="R76", tries=p, den_ruling="R76")
    assert rec["name"] == "T3_c2_s2_s05" and rec["den_ruling"]["row"] == "R76"


@pytest.mark.parametrize("g,s", [
    (None, None), (None, 0.85), (None, 0.5), (None, 0.67), (1.5, None), (0.5, None), (2.0, 0.67), (INF, None),
    (INF, 0.67), (2.0, None),  # row R76's, R47's and line 3's points: not row R79's
    (1.0, 0.35), (None, 0.350001), (None, 0.3499999), (None, 0.36), (None, 0.3), (0.5, 0.35), (INF, 0.35),
    (0.5, None), (0.5, 0.500001), (0.500001, 0.5), (2.000001, 0.5), (1.5, 0.67), (2.0, 0.85), (1.0, 0.5),
    (-INF, 0.5), (float("nan"), 0.5), (INF, float("nan")), (INF, INF), (INF, -0.5), (3.0, 0.5), (0.0, 0.5),
    (None, 1.0), (1.0, 1.0), (True, 0.5), (INF, True), ("x", 0.5), ("in f", 0.5), (None, "0,35"), (None, ""),
])
def test_a_near_miss_value_is_refused(tmp_path, g, s):
    """Every command near one of the five (a value off, an option added, dropped or swapped, row R76's, row R47's and
    line 3's own points) is refused, whatever --den-ruling says."""
    p = ledger(tmp_path)
    for den in (None, "R79"):
        with pytest.raises(bk.Refused, match="(not one of row R79's five option strings|is refused under)"):
            bk.check_r79("c2-s2", g, s, "R79", p, den_ruling=den)


def test_value_equal_spellings_are_the_same_value_as_on_every_line(tmp_path):
    """The command line reads numbers (typer's float), as every earlier line: 0.350 is the value 0.35 and Inf the value
    inf, so they name the same file; the row and the launcher match the strings themselves (tested there)."""
    p = ledger(tmp_path)
    assert bk.check_r79("c2-s2", None, "0.350", "R79", p, den_ruling="R79")["name"] == "T3_c2_s2_s035"
    assert bk.check_r79("c2-s2", "Inf", "0.50", "R79", p, den_ruling="R79")["name"] == "T3_c2_s2_gonly_s05"


@pytest.mark.parametrize("form", [f for f in bk.FORMS if f != "c2-s2"] + ["", None, "c2-s2 "])
def test_another_form_is_refused(tmp_path, form):
    for name in NAMES:
        with pytest.raises(bk.Refused, match="on --form c2-s2 only"):
            bk.check_r79(form, captain_ruling="R79", tries=ledger(tmp_path), den_ruling="R79", **KW[name])


@pytest.mark.parametrize("extra", [{"composition_from": "allowed-s1"}, {"depth_from": "allowed-s1"},
                                   {"depth_exponent": 0.5}, {"depth_form": "global"}])
def test_the_composition_and_depth_options_beside_the_ruling_are_refused(tmp_path, extra):
    for name in NAMES:
        with pytest.raises(bk.Refused, match="--captain-ruling R79 takes no --(composition-from|depth-)"):
            bk.check_r79("c2-s2", captain_ruling="R79", tries=ledger(tmp_path), den_ruling="R79", **KW[name], **extra)


@pytest.mark.parametrize("den", [None, "R76", "R47", "R22", "R30", "r79", "R79 ", "R790", "", "R7"])
@pytest.mark.parametrize("name", NAMES)
def test_each_string_is_taken_only_with_den_ruling_r79(tmp_path, name, den):
    with pytest.raises(bk.Refused, match="is taken only with --den-ruling R79"):
        bk.check_r79("c2-s2", captain_ruling="R79", tries=ledger(tmp_path), den_ruling=den, **KW[name])


def test_a_string_the_table_did_not_lift_would_take_no_den_ruling(tmp_path, monkeypatch):
    """The table's lift column is read as row R76's is: were a string not lifted, --den-ruling would be refused for it
    (no row R79 string is such: every one carries a size scale)."""
    opt, g, s, _ = bk.R79_FILES["T3_c2_s2_s035"]
    monkeypatch.setitem(bk.R79_FILES, "T3_c2_s2_s035", (opt, g, s, False))
    with pytest.raises(bk.Refused, match="does not lift the 5% check for this string"):
        bk.check_r79("c2-s2", captain_ruling="R79", tries=ledger(tmp_path), den_ruling="R79", size_scale=0.35)


def ok_row_but(**kw) -> str:
    return ROW_R22 + ROW_R47 + ROW_R76 + r79_row(**kw)


ROW_REFUSALS = [
    (ROW_R22 + ROW_R47 + ROW_R76, "no row R79"),
    (ROWS + ROW_R79, "two rows R79"),
    (ok_row_but().replace("| R79 (", "| R79b ("), "another row id"),
    (ok_row_but(opener="R79 (the captain's ruling of a day"), "not the captain's order"),
    (ok_row_but(opener="R79 (the captain's named exception of a day"), "not the captain's order"),
    (ok_row_but(ruling="a run"), "no ruling, no run"),
    (ok_row_but(captain="--captain-ruling R79"), "the ruling without backticks"),
    (ok_row_but(captain="`--captain-ruling R76`"), "another ruling"),
    (ok_row_but(den="--den-ruling R79"), "the lift's ruling without backticks"),
    (ok_row_but(den="`--den-ruling R76`"), "another lifting row"),
    (ok_row_but(words="the 5% check gates nothing for them"), "the lifting words changed"),
    (ok_row_but(words="the 5% check gates nothing"), "the lifting words cut"),
    (ok_row_but(quote=QUOTE.replace("gates nothing.\"", "gates nothing for it.\"")), "the lift-phrase trap"),
    (ok_row_but(lifted=NAMES[:4]), "the lifting part misses one string"),
    (ok_row_but(lifted=NAMES[1:]), "the lifting part misses the first string"),
    (ROWS.replace(" (`--den-ruling R79`); decision", " (`--den-ruling R79`); the 5% check gates nothing for it; "
                  "decision"), "two parts say the words"),
    (ok_row_but(names=NAMES[:4]), "a string missing from the table"),
    (ok_row_but(pairs=", ".join(f"`{OPT[n]}` ({m})" for n, m in zip(NAMES, [*NAMES[1:], NAMES[0]], strict=True))),
     "the strings paired with other files"),
    (ok_row_but(pairs=", ".join(f"`{OPT[n]}`" for n in NAMES)), "the strings without their files"),
    (ok_row_but(pairs=", ".join(f"{OPT[n]} ({n})" for n in NAMES)), "the strings without backticks"),
    (ok_row_but(pairs=", ".join(f"`a`{OPT[n]}`b` ({n})" for n in NAMES)), "not as whole strings"),
    (ok_row_but(pairs=", ".join(f"`{OPT[n]} ` ({n})" for n in NAMES)), "a trailing space inside the backticks"),
    (ok_row_but(pairs=", ".join(f"`{OPT[n].replace('c2-s2', 'c2-s1')}` ({n})" for n in NAMES)), "another form"),
    # the near-miss strings: a row that spells a string otherwise opens no file
    (ROWS.replace("--size-scale 0.35`", "--size-scale 0.350`"), "0.35 written 0.350"),
    (ROWS.replace("--gata4-exponent inf --size-scale 0.5`", "--gata4-exponent Inf --size-scale 0.5`"),
     "inf written Inf"),
    (ROWS.replace("--gata4-exponent 2.0 --size-scale 0.5`", "--gata4-exponent 2 --size-scale 0.5`"), "2.0 written 2"),
    ("| R79 (the captain's order of a day, ruling, no run; not a try) | T3 | "
     + ", ".join(f"`{OPT[n]}` ({n})" for n in NAMES) + " `--captain-ruling R79`; for each of "
     + ", ".join(f"`{OPT[n]}`" for n in NAMES) + ", the 5% check gates nothing for it (`--den-ruling R79`) |\n",
     "the table outside the first cell"),
]


@pytest.mark.parametrize("rows,why", ROW_REFUSALS, ids=[w for _, w in ROW_REFUSALS])
def test_a_row_r79_that_does_not_hold_the_table_opens_no_file(tmp_path, rows, why):
    p = ledger(tmp_path, rows)
    for name in NAMES:
        with pytest.raises(bk.Refused, match="has no single row R79 whose first cell is the captain's order"):
            bk.check_r79("c2-s2", captain_ruling="R79", tries=p, den_ruling="R79", **KW[name])


def test_the_row_check_is_build_outside_form_s_plus_the_table(tmp_path, monkeypatch):
    bo = bk._outside_form()
    seen = []
    real = bo.captain_order_names
    monkeypatch.setattr(bo, "captain_order_names", lambda row, option, path=None: seen.append((row, option, path))
                        or real(row, option, path))
    p = ledger(tmp_path)
    for name in NAMES:
        seen.clear()
        assert bk.check_r79("c2-s2", captain_ruling="R79", tries=p, den_ruling="R79", **KW[name])["name"] == name
        assert seen == [("R79", f"`{OPT[name]}`", p)]


@pytest.mark.parametrize("ruling", [None, "R22", "R47", "R41", "R53", "R36", "R30", "R76", "r79", "R79 "])
@pytest.mark.parametrize("name", NAMES)
def test_each_string_under_any_other_ruling_is_refused_as_before(world, tmp_path, monkeypatch, ruling, name):
    """Under --captain-ruling R22 only 0.67 is a size scale; under R76 only row R76's strings; under any other ruling
    every line refuses; R79's own code is never reached."""
    for fn in LINE_FUNCS:
        monkeypatch.setattr(bk, fn, lambda *a, _n=fn, **k: pytest.fail(f"{_n} reached under {ruling!r}"))
    out = tmp_path / "x.h5ad"
    for den in (None, "R79"):
        with pytest.raises(bk.Refused):
            build(world, out, captain_ruling=ruling, tries=ledger(tmp_path), den_ruling=den, **KW[name])
        assert nothing_written(out) and loader_restored()
    if ruling == "R76":
        with pytest.raises(bk.Refused, match="not one of row R76's five option strings"):
            build(world, out, captain_ruling="R76", tries=ledger(tmp_path), **KW[name])
    if ruling == "R22":
        with pytest.raises(bk.Refused, match="only 0.67"):  # every R79 string's size scale is another
            build(world, out, captain_ruling="R22", tries=ledger(tmp_path), **KW[name])


@pytest.mark.parametrize("name", list(R76_OPT))
def test_each_r76_string_under_r79_is_refused_and_builds_as_before_under_r76(world, runs, tmp_path, monkeypatch, name):
    for den in (None, "R79", "R76"):
        out = tmp_path / "x.h5ad"
        with pytest.raises(bk.Refused, match="not one of row R79's five option strings"):
            build(world, out, captain_ruling="R79", tries=ledger(tmp_path), den_ruling=den, **R76_KW[name])
        assert nothing_written(out)
    r = use(runs, monkeypatch, "mild")
    lift = {"den_ruling": "R76"} if name in R76_LIFTED else {}
    side = build(world, tmp_path / f"{name}.h5ad", run_dir=r.run, captain_ruling="R76", tries=ledger(tmp_path),
                 **R76_KW[name], **lift)
    assert side["combined_build"]["r76_line"]["file"] == f"{name}.h5ad" and "r79_line" not in side["combined_build"]


def test_den_ruling_r79_lifts_no_f_file_and_no_r76_file(world, tmp_path):
    """Row R79 names none of F1, F2, F3's whole strings: --den-ruling R79 is refused for them (check_den_lift), and
    row R76 takes only --den-ruling R76."""
    p = ledger(tmp_path)
    for kw in ({"gata4_exponent": INF}, {"size_scale": 0.67}, {"gata4_exponent": INF, "size_scale": 0.67}):
        with pytest.raises(bk.Refused, match="--den-ruling R79 is refused"):
            bk.check_c2_gonly("c2-s2", captain_ruling="R22", tries=p, den_ruling="R79", **kw)
    for name in R76_OPT:
        with pytest.raises(bk.Refused, match="(is taken only with --den-ruling R76|is refused for)"):
            bk.check_r76("c2-s2", captain_ruling="R76", tries=p, den_ruling="R79", **R76_KW[name])


# ------------------------------------------------------------------------------------------- the files at the build


@pytest.mark.parametrize("name", NAMES)
def test_each_file_is_c2_s_command_with_its_factor_and_records_both_checks(world, runs, tmp_path, monkeypatch, name):
    r = use(runs, monkeypatch, "mild")
    spy = Spy(monkeypatch)
    c2 = build(world, tmp_path / "c2.h5ad", run_dir=r.run)
    out = tmp_path / f"{name}.h5ad"
    side = build_r(world, tmp_path, name, out, run_dir=r.run)
    assert len(spy.seen) == 2
    g, s = POINTS[name]
    want, rec = bk.r79_factor(r.m, r.g, g, s)
    assert np.array_equal(spy.seen[1], want) and np.array_equal(spy.seen[0], bk.combine("c2", 1.0, r.m, r.g)[0])
    # C2's command with its own option string: spec, inputs, streams, composition, seed, count and board are C2's
    assert side["spec"] == c2["spec"] and side["inputs"] == c2["inputs"] and side["spec_source"] == c2["spec_source"]
    assert side["inputs"]["rng_streams"] == ["submission:barred:s2", "submission:barred:s2:assemble"]
    assert side["construction"]["composition"] == c2["construction"]["composition"]
    assert side["construction"]["factor_exponent"] == 1.0 and side["output"]["check"]["ok"]
    assert side["construction"]["depth_ratios"] is None and set(side["spec"]) == {"factor", "composition"}
    cb, cb2 = side["combined_build"], c2["combined_build"]
    assert {k: v for k, v in cb.items() if k not in ("combination", "r79_line")} == {
        k: v for k, v in cb2.items() if k != "combination"}
    assert cb["form"] == "c2-s2" and cb["combination"] == rec and cb["depth_step"] is None
    assert "r76_line" not in cb and "c2_gonly_line" not in cb
    line = cb["r79_line"]
    assert line["captain_ruling"]["row"] == "R79" and line["captain_ruling"]["name"] == name
    assert line["file"] == f"{name}.h5ad" and line["option"] == f"`{OPT[name]}`" and line["line"] == bk.R79_LINE
    assert line["g"] == ("inf" if g == INF else g) and line["size_scale"] == s and line["c"] == rec["c"]
    assert line["captain_ruling"]["g"] == line["g"]
    assert line["powers_applied"] == {bk.TABLE_M: rec["mab21l2_exponent"], bk.TABLE_G: rec["gata4_weight"]}
    assert line["closed_lever"] == bk.R79_CLOSED_LEVER and line["list_entry"] == bk.R79_LIST_ENTRY
    assert line["den_rule"] == bk.R79_DEN_RULE and line["functions"] == bk.r79_functions(g)["branch"]
    assert line["c2_reference_sha256"] == "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"
    sha = hashlib.sha256(np.ascontiguousarray(want).tobytes()).hexdigest()
    assert line["factor_sha256_applied"] == line["guards"]["factor_sha256"] == sha
    gd = line["guards"]
    assert "genes" not in gd and gd["n_genes"] == G and gd["tables_sha256"] == r.pins
    assert gd["passed"] and gd["range"]["passed"] and gd["den_check"]["passed"] and gd["den_check"]["within"]
    assert line["den_ruling"]["row"] == "R79" and gd["den_check"]["lifted_by"] == line["den_ruling"]
    assert gd["den_check"]["limit"]["status"] == "lifted by row R79: the 5% check gates nothing for it"
    want_g = bk.r79_functions(g)["guards"](r.m, r.g, g, s, den_of(world), line["den_ruling"])
    assert as_json({k: gd[k] for k in want_g}) == as_json(want_g)
    assert np.array_equal(np.array(gd["den"]), den_of(world)) and gd["wt_sha256"] == t3d.sha256(world.wt)
    on_disk = json.loads(out.with_suffix(".json").read_text())["combined_build"]["r79_line"]
    assert on_disk["file"] == f"{name}.h5ad" and on_disk["guards"]["range"] == gd["range"]
    assert on_disk["g"] == ("inf" if g == INF else g)  # a string in the record, never a bare non-finite number
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


def test_the_five_files_are_five_factors_none_c2_s_or_an_earlier_file_s(world, runs, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    shas = {name: bk.factor_sha256(bk.r79_factor(r.m, r.g, *POINTS[name])[0]) for name in NAMES}
    earlier = {bk.factor_sha256(bk.combine("c2", 1.0, r.m, r.g)[0])}
    earlier |= {bk.factor_sha256(bk.r76_factor(r.m, r.g, 1.0 if g is None else g, 1.0 if s is None else s)[0])
                for _, g, s, _ in bk.R76_FILES.values()}
    earlier |= {bk.factor_sha256(bk.c2_gonly_factor(r.m, r.g, *p)[0]) for p in bk.GONLY_POINTS.values()}
    earlier.add(bk.factor_sha256(bk.c2_gata4_weighted(r.m, r.g, 2.0)[0]))
    assert len(set(shas.values())) == 5 and not set(shas.values()) & earlier


@pytest.mark.parametrize("name", NAMES)
def test_the_5_percent_check_gates_nothing_for_each(world, runs, tmp_path, monkeypatch, name):
    """FAR: every file's den-weighted total is far from C2's and its range within C2's bounds; each builds and records
    the figures, which gate nothing."""
    r = use(runs, monkeypatch, "far")
    g, s = POINTS[name]
    pre = bk.r79_functions(g)["guards"](r.m, r.g, g, s, den_of(world))
    assert not pre["den_check"]["within"] and pre["range"]["passed"] and not pre["passed"]  # the premise
    side = build_r(world, tmp_path, name, tmp_path / f"{name}.h5ad", run_dir=r.run)
    d = side["combined_build"]["r79_line"]["guards"]["den_check"]
    assert d["passed"] and not d["within"] and abs(d["den_shift_vs_c2"]) > 0.05
    assert d["lifted_by"]["row"] == "R79" and side["output"]["check"]["ok"]


@pytest.mark.parametrize("name", ["T3_c2_s2_s035", "T3_c2_s2_gonly_s05"])
def test_the_range_guard_stops_a_file_before_any_cell_is_drawn(world, runs, tmp_path, monkeypatch, name):
    r = use(runs, monkeypatch, "mild")
    c2, _ = bk.combine("c2", 1.0, r.m, r.g)
    wide = np.ones(G)
    wide[0] = 3.0 * c2.max()
    fn = "c2_gonly_factor" if POINTS[name][0] == INF else "r76_factor"
    monkeypatch.setattr(bk, fn, lambda *a, **k: (wide, {"g": 2.0, "size_scale": 0.5, "c": 1.0,
                                                        "mab21l2_exponent": 0.5, "gata4_weight": 1.0}))
    monkeypatch.setattr(bk.bs, "construct", lambda *a, **k: pytest.fail("a cell was drawn past a failed check"))
    out = tmp_path / f"{name}.h5ad"
    with pytest.raises(bk.R79GuardFailed, match="the range guard"):
        build_r(world, tmp_path, name, out, run_dir=r.run)
    assert nothing_written(out) and loader_restored()
    assert issubclass(bk.R79GuardFailed, bk.GonlyGuardFailed)  # the command line's exit 6


@pytest.mark.parametrize("name", ["T3_c2_s2_g15_s05", "T3_c2_s2_gonly_s05"])
def test_a_factor_other_than_the_checks_exits_5_and_writes_nothing(world, runs, tmp_path, monkeypatch, name):
    r = use(runs, monkeypatch, "mild")
    fn = "run_gonly_guards" if POINTS[name][0] == INF else "run_r76_guards"
    real = getattr(bk, fn)
    monkeypatch.setattr(bk, fn, lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.FactorMismatch, match="not the one the guards describe"):
        build_r(world, tmp_path, name, out, run_dir=r.run)
    assert nothing_written(out)
    monkeypatch.setattr(bk, fn, lambda *a, **k: real(*a, **k) | {"genes": ["Z"]})
    with pytest.raises(bk.FactorMismatch, match="genes"):
        build_r(world, tmp_path, name, out, run_dir=r.run)
    assert nothing_written(out) and loader_restored()


def test_the_checks_read_the_pinned_tables_only(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    monkeypatch.setattr(bk, "TABLE_SHA256", {**r.pins, bk.TABLE_G: "0" * 64})
    for name in ("T3_c2_s2_g05_s05", "T3_c2_s2_gonly_s05"):
        out = tmp_path / f"{name}.h5ad"
        with pytest.raises(bk.Refused, match="pinned table"):
            build_r(world, tmp_path, name, out, run_dir=r.run)
        assert nothing_written(out)


def test_an_existing_output_is_refused_before_the_checks(world, runs, tmp_path, monkeypatch):
    r = use(runs, monkeypatch, "mild")
    out = tmp_path / "T3_c2_s2_s035.h5ad"
    out.write_bytes(b"x")
    for fn in ("run_r76_guards", "run_gonly_guards"):
        monkeypatch.setattr(bk, fn, lambda *a, **k: pytest.fail("the checks ran past an existing file"))
    with pytest.raises(bk.Refused, match="exists; pass --overwrite"):
        build_r(world, tmp_path, "T3_c2_s2_s035", out, run_dir=r.run)
    assert out.read_bytes() == b"x"


def test_an_unclean_tree_is_refused(world, tmp_path, own_repo):
    (own_repo / "stray.py").write_text("x = 1\n")
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match="not clean"):
        build_r(world, tmp_path, "T3_c2_s2_g2_s05", out)
    assert nothing_written(out)


# ------------------------------------------------------------------------------- every earlier command as it was


OLD_SETS = [(f, {}, "world") for f in bk.FORMS] + [
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "world"),
    ("c2-s2", {"size_scale": 0.67, "captain_ruling": "R22", "den_ruling": "R47"}, "mild"),
    ("c2-s2", {"gata4_exponent": INF, "captain_ruling": "R22"}, "mild"),
    ("c2-s2", {"gata4_exponent": INF, "size_scale": 0.67, "captain_ruling": "R22", "den_ruling": "R47"}, "mild"),
    *[("c2-s2", {**R76_KW[n], "captain_ruling": "R76", **({"den_ruling": "R76"} if n in R76_LIFTED else {})}, "mild")
      for n in R76_OPT],
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
    assert plain["combined_build"] == trapped["combined_build"] and "r79_line" not in plain["combined_build"]


# --------------------------------------------------------------------------------------------------------- CLI


def test_all_five_through_the_typer_command_line(world, runs, tmp_path, monkeypatch):
    """Each of the five by its own command line, as the launcher passes it: exit 0, the R79 summary line (g = inf for
    T3_c2_s2_gonly_s05, the record's string printed as it is), the file and its sidecar."""
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    boards = th.cli_boards(tmp_path)
    r = use(runs, monkeypatch, "mild")
    runner = CliRunner()
    for name in NAMES:
        ok = tmp_path / f"{name}.h5ad"
        res = runner.invoke(bk.app, ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic",
                                     "--boards-dir", str(boards), "--run-dir", str(r.run), "--seed", str(SEED),
                                     "--threads", "2", *OPT[name].split(), "--captain-ruling", "R79",
                                     "--den-ruling", "R79", "--out", str(ok)])
        assert res.exit_code == 0, (name, res.output)
        g, s = POINTS[name]
        gtxt = "inf" if g == INF else f"{g:g}"
        assert (f"row R79's refinement round 4 on C2's command: {name}.h5ad (g = {gtxt}, s = {s:g}); the range guard "
                "passed; the 5% den check gates nothing for it (row R79)") in res.output
        assert res.output.splitlines()[-1] == f"OK: {ok}" and "row R76's" not in res.output
        line = json.loads(ok.with_suffix(".json").read_text())["combined_build"]["r79_line"]
        assert line["file"] == f"{name}.h5ad" and line["captain_ruling"]["name"] == name
        assert line["g"] == ("inf" if g == INF else g)
    assert bk.r79_g_text("inf") == "inf" and bk.r79_g_text(2.0) == "2" and bk.r79_g_text(0.5) == "0.5"


def test_the_cli_refusals_and_stops(world, runs, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    boards = th.cli_boards(tmp_path)
    r = use(runs, monkeypatch, "mild")

    def args(run, *extra, out):
        return ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
                str(boards), "--run-dir", str(run), "--seed", str(SEED), "--threads", "2", *extra, "--out", str(out)]

    runner = CliRunner()
    refusals = [[*OPT[n].split(), "--captain-ruling", "R22"] for n in NAMES]  # each string under R22: as before
    refusals += [[*OPT[n].split(), "--captain-ruling", "R76", "--den-ruling", "R76"] for n in NAMES]
    refusals += [[*OPT[n].split(), "--captain-ruling", "R79"] for n in NAMES]  # no --den-ruling R79
    refusals += [[*OPT[n].split(), "--captain-ruling", "R79", "--den-ruling", "R76"] for n in NAMES]
    refusals += [[*OPT[n].split(), "--den-ruling", "R79"] for n in NAMES]  # no ruling
    refusals += [[*R76_OPT[n].split(), "--captain-ruling", "R79", "--den-ruling", "R79"] for n in R76_OPT]
    refusals += [["--form", "c2-s2", "--captain-ruling", "R79", "--den-ruling", "R79"],
                 ["--form", "c2-s2", "--size-scale", "0.351", "--captain-ruling", "R79", "--den-ruling", "R79"],
                 ["--form", "c2-s2", "--gata4-exponent=-inf", "--size-scale", "0.5", "--captain-ruling", "R79",
                  "--den-ruling", "R79"],
                 ["--form", "c2-s2", "--gata4-exponent", "inf", "--size-scale", "0.67", "--captain-ruling", "R79",
                  "--den-ruling", "R79"],
                 ["--form", "c2-s2", "--gata4-exponent", "1.0", "--size-scale", "0.35", "--captain-ruling", "R79",
                  "--den-ruling", "R79"],
                 ["--form", "c1-s2", "--gata4-exponent", "1.5", "--size-scale", "0.5", "--captain-ruling", "R79",
                  "--den-ruling", "R79"],
                 ["--form", "c2-s2", "--size-scale", "0.35", "--captain-ruling", "R79", "--den-ruling", "R79",
                  "--depth-from", "allowed-s1"],
                 ["--form", "c2-s2", "--size-scale", "0.35", "--captain-ruling", "R79", "--den-ruling", "R79",
                  "--composition-from", "allowed-s1"]]
    for extra in refusals:
        out = tmp_path / "r.h5ad"
        res = runner.invoke(bk.app, args(r.run, *extra, out=out))
        assert res.exit_code == 2 and "refused" in res.output, (extra, res.output)
        assert nothing_written(out)
    # the range guard stops the file: exit 6, nothing written
    c2, _ = bk.combine("c2", 1.0, r.m, r.g)
    wide = np.ones(G)
    wide[0] = 3.0 * c2.max()
    real_gonly = bk.c2_gonly_factor
    monkeypatch.setattr(bk, "c2_gonly_factor", lambda *a, **k: (wide, {"g": "inf", "size_scale": 0.5, "c": 1.0,
                                                                       "mab21l2_exponent": 0.0, "gata4_weight": 0.5}))
    bad = tmp_path / "stop.h5ad"
    res = runner.invoke(bk.app, args(r.run, *OPT["T3_c2_s2_gonly_s05"].split(), "--captain-ruling", "R79",
                                     "--den-ruling", "R79", out=bad))
    assert res.exit_code == 6 and "refused (the range guard or the 5% den check" in res.output, res.output
    assert nothing_written(bad)
    monkeypatch.setattr(bk, "c2_gonly_factor", real_gonly)
    # the factor is not the checks': exit 5, nothing written
    real = bk.run_r76_guards
    monkeypatch.setattr(bk, "run_r76_guards", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    mis = tmp_path / "mismatch.h5ad"
    res = runner.invoke(bk.app, args(r.run, *OPT["T3_c2_s2_g2_s05"].split(), "--captain-ruling", "R79",
                                     "--den-ruling", "R79", out=mis))
    assert res.exit_code == 5 and nothing_written(mis)
    monkeypatch.setattr(bk, "run_r76_guards", real)
    # a ledger without row R79: refused
    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path, ROW_R22 + ROW_R47 + ROW_R76, name="no_r79.md"))
    res = runner.invoke(bk.app, args(r.run, *OPT["T3_c2_s2_s035"].split(), "--captain-ruling", "R79",
                                     "--den-ruling", "R79", out=tmp_path / "n.h5ad"))
    assert res.exit_code == 2 and "no single row R79" in res.output and nothing_written(tmp_path / "n.h5ad")


def test_the_r79_check_comes_first_on_the_command_line(world, runs, tmp_path, monkeypatch):
    """With --captain-ruling R79 no other line's check runs before R79's (the options' line would refuse a size
    scale other than 0.67, row R76's would refuse every R79 string)."""
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    r = use(runs, monkeypatch, "mild")
    order = []
    for fn in ("check_r79", "check_r76", "check_depth_line", "check_c2_gonly", "check_composition_from",
               "check_c2_gata4", "check_gata4"):
        real = getattr(bk, fn)
        monkeypatch.setattr(bk, fn, lambda *a, _f=fn, _r=real, **k: order.append(_f) or _r(*a, **k))
    res = CliRunner().invoke(bk.app, ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic",
                                      "--boards-dir", str(th.cli_boards(tmp_path)), "--run-dir", str(r.run),
                                      "--threads", "2", *OPT["T3_c2_s2_s035"].split(), "--captain-ruling", "R79",
                                      "--den-ruling", "R79", "--out", str(tmp_path / "f.h5ad")])
    assert res.exit_code == 0, res.output
    assert order == ["check_r79", "check_r79"]  # the command line's, then build_form's: no other line's
