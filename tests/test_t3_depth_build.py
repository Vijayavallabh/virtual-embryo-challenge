"""The depth line in experiments/t3_decomp/build_two_knockout.py (the T3 ERROR design of the design round of 2026-11-25,
section 3, proposal P1, dD, as its cold challenge fixes it; row R53, the captain's order of 2026-11-24): on
`--form c2-s2`, `--depth-from allowed-s1` (dd, T3_c2_s2_dd.h5ad), with `--depth-exponent 0.5` (ddh, T3_c2_s2_ddh.h5ad)
or `--depth-form global` (dg, T3_c2_s2_dg.h5ad), each only with `--captain-ruling R53` while row R53's first cell is the
captain's order naming the file's option string and `--form c2-s2`, each as one whole backtick-delimited string. The
frozen allowed branch's depth step is added to the barred s2's spec where the frozen builder resolves it, and its
per-label ratios (build_submission.depth_ratios on the pinned decision file's submissions.parameters) enter as
construct's `ratio`, after C2's factor on the same generator.

Held here: the ratios by hand (as they are, to the power 0.5, one cell-weighted geometric mean over all labels, 1.0
where none); the refusals; the ruling-row check (build_outside_form.captain_order_names plus the whole-string matches),
on synthetic ledgers and on row R53's first cell as the repository's ledger holds it, read by code; without the options
every byte of C2 and of every other command as it was, the line's code never reached; with them C2's cells, rows,
factor draws, composition, streams and coordinates, only the expression changed; dd byte for byte the frozen builder's
own depth step on C2's spec; the counts the depth step receives equal to C2's final counts (across two builds, and
inside the build, where the check stops the file, exit 7, before anything is written); the cells whose thinning draw
would be empty kept unthinned as the frozen thin keeps them, counted on the very draw and recorded; two builds, one
sha256; the command line.

Offline, synthetic and fast: the throwaway world of tests/test_t3_two_knockout_build.py (a small synthetic WT file,
synthetic `transfer_mab21l2_all` and `gata4_pooled` tables), a synthetic decision file with the allowed s1 entry and
its depth ratios, and small ledger files under tmp_path. Nothing here builds from released data or pins bytes built on
one machine; of the repository's ledger only row R53's first cell is read (by code, never printed), and of the
committed decision file only its sha256 and the shape of the entry the line reads.
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
COMMITTED_DECISION = ROOT / "experiments" / "t3_decomp" / "results" / "gata4" / "gata4_decisions.json"
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
L2 = 1.3826521560990648
BLOCK = json.loads(json.dumps(bs.composition_spec("central", L2, tx.RULES, "Gata4 replicates (D7)")))
# the synthetic frozen ratios: below and above 1 (thinning and top-up), and one label (Endo) left at 1.0 as the frozen
# rule leaves a label with too few cells
POOLED = {"vCM1": 0.6, "OFT/RV-CM": 0.75, "aSHF": 0.5, "pSHF": 0.45, "aPHM": 0.46, "pPHM": 0.55, "NCC": 0.56,
          "SE": 2.0, "Foregut": 0.62, "Peri": 1.4}
TINY = {**POOLED, "NCC": 0.004}  # NCC thinned to almost nothing: about half its cells' draws come out empty
ARMS = ("dd", "ddh", "dg")
OPT = {"dd": "--depth-from allowed-s1", "ddh": "--depth-from allowed-s1 --depth-exponent 0.5",
       "dg": "--depth-from allowed-s1 --depth-form global"}
FORM = "--form c2-s2"
FILES = {"dd": "T3_c2_s2_dd.h5ad", "ddh": "T3_c2_s2_ddh.h5ad", "dg": "T3_c2_s2_dg.h5ad"}
KW = {"dd": {"depth_from": "allowed-s1"}, "ddh": {"depth_from": "allowed-s1", "depth_exponent": 0.5},
      "dg": {"depth_from": "allowed-s1", "depth_form": "global"}}
LINE_FUNCS = ("depth_parameters", "depth_ratios_applied", "depth_step", "depth_checks", "_twin_generator")
OTHER_ROWS = ("| R22 (the captain's order of a day; ruling, no run, not a try; the builder is to take "
              "`--gata4-exponent 2.0` on `--form c1-s2`, `--form c2-s2 --gata4-exponent 2.0`, "
              "`--form c2-s2 --gata4-exponent inf`, `--form c2-s2 --size-scale 0.67` and "
              "`--form c2-s2 --gata4-exponent inf --size-scale 0.67`, each only with `--captain-ruling R22`) "
              "| T3 | x |\n"
              "| R36 (the captain's named exception of a day; ruling, no run, not a try; the builder is to take "
              "`--composition-from allowed-s1` with `--form c2-s2` only with `--captain-ruling R36`) | T3 | y |\n")


def row(*arms: str, rid: str = "R53", opener: str = "the captain's order of a day, ruling, no run",
        form: str = f"`{FORM}`", ruling: str | None = None, note: bool | None = None) -> str:
    """A ruling row in row R53's form: the files' option strings, the form and the ruling, each in backticks. R53's
    note on the two variants names `--depth-from allowed-s1` itself, so it is kept only for a row naming all three
    files (note=None) unless asked for."""
    names = ", ".join(f"`{OPT[a]}`" for a in arms)
    ruling = f"`--captain-ruling {rid}`" if ruling is None else ruling
    note = set(arms) == set(ARMS) if note is None else note
    aside = (" (the options `--depth-exponent 0.5` and `--depth-form global` only beside `--depth-from allowed-s1`)"
             if note else "")
    return (f"| {rid} ({opener}: \"Yes: my named exception on T3's depth\"; the T3 two-knockout builder is to take "
            f"{names}{aside} on {form} only with {ruling}; not a try) | T3 | z |\n")


ROW = row(*ARMS)


def ledger(tmp_path: Path, rows: str = ROW, name: str = "TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text(f"# a throwaway ledger\n| R21 (another row) | T2 | y |\n{OTHER_ROWS}{rows}| R41 (a yes) | T3 | z |\n")
    return p


def write_decision(path: Path, pooled: dict | None = POOLED, depth=True, parameters: dict | None = None) -> Path:
    """A synthetic decision file: D7m (the barred branch), D7 (row R36's m), the allowed s1 entry (its depth flag as
    given; none when None) and submissions.parameters with the depth ratios (or `parameters` as given)."""
    rules = {"D7m": {"outcome": True, "m": tx.M7M, "n_counted": 1}, "D7": {"outcome": True, "m": L2}}
    ent = {"factor": {"source": "gata4", "exponent": 1.0}, "composition": BLOCK, "stage": False, "coverage": False}
    if depth is not None:
        ent["depth"] = depth
    params = parameters if parameters is not None else ({} if pooled is None else {"depth_ratio_pooled": pooled})
    return tx._decision(path, rules=rules, submissions={"allowed": {"s1": ent}, "parameters": params})


@pytest.fixture(scope="module")
def dw(world):
    """The line's decision files beside the world's (which carries D7m only), and the reference's vocabulary as the
    build loads it (load_reference on the synthetic WT file: the labels and WT cells per label the ratios use)."""
    dec = write_decision(world.d / "gata4_decisions_depth.json")
    tiny = write_decision(world.d / "gata4_decisions_depth_tiny.json", TINY)
    ref = bs.load_reference(world.wt, "gata4", tx.RULES, None, 2)
    return SimpleNamespace(dec=dec, sha=t3d.sha256(dec), tiny=tiny, tiny_sha=t3d.sha256(tiny), V=ref.V)


@pytest.fixture(autouse=True)
def pins(dw, monkeypatch):
    """The synthetic decision file pinned for the line and for row R36 (the committed pin is tested below)."""
    monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", dw.sha)
    monkeypatch.setattr(bk, "R36_DECISION_SHA256", dw.sha)


def build(world, dec: Path, out: Path, form: str = "c2-s2", **kw) -> dict:
    return bk.build_form(dec, form, world.wt, kw.pop("board", world.board), N, out,
                         run_dir=kw.pop("run_dir", world.run), seed=SEED, threads=2, **kw)


def build_arm(world, dw, tmp_path: Path, arm: str, out: Path, rows: str = ROW, dec: Path | None = None, **kw) -> dict:
    return build(world, dec or dw.dec, out, captain_ruling=kw.pop("captain_ruling", "R53"),
                 tries=ledger(tmp_path, rows, name=f"{out.stem}_TRIES.md"), **KW[arm], **kw)


def nothing_written(out: Path) -> bool:
    return not out.exists() and not out.with_suffix(".json").exists() and not Path(f"{out}.stamp.json").exists()


def hooks_restored() -> bool:
    """The builder's hooks are gone (a test's own spies, if any, may stay)."""
    frozen = {"resolve_spec": bs.resolve_spec, "depth_ratios": bs.depth_ratios, "construct": bs.construct,
              "load_factor": bs.load_factor, "inputs_for": bs.inputs_for}
    return (all(f.__module__ == "build_submission" and f.__name__ == n for n, f in frozen.items())
            and not t3d.lognorm.__qualname__.startswith("_lognorm_spy")
            and not t3d.rescale_depth.__qualname__.startswith("_rescale_spy"))


def by_hand(V, pooled: dict, arm: str) -> np.ndarray:
    """The ratios the arm applies, written out: r_t from the pooled table (1.0 for a label without one); dd r_t; ddh
    sqrt(r_t); dg one value exp(sum_t n_t ln r_t / sum_t n_t) over all labels, n_t the WT cells of label t."""
    r = np.array([pooled.get(t, 1.0) for t in V.labels], dtype=np.float64)
    if arm == "dd":
        return r
    if arm == "ddh":
        return np.sqrt(r)
    n = np.asarray(V.count, dtype=np.float64)
    return np.full(len(r), math.exp(float(np.sum(n * np.log(r)) / np.sum(n))))


def read(path: Path):
    import anndata as ad

    return ad.read_h5ad(path)


class Observed:
    """What the build's t3_decompose.build hands to rescale_depth and to lognorm (calls inside build only)."""

    def __init__(self, monkeypatch):
        self.rescale: list = []
        self.lognorm: list = []
        self.inside = False
        real_build, real_rescale, real_lognorm = t3d.build, t3d.rescale_depth, t3d.lognorm

        def spy_build(*a, **k):
            self.inside = True
            try:
                return real_build(*a, **k)
            finally:
                self.inside = False

        def spy_rescale(C, ratio, P, rng):
            res = real_rescale(C, ratio, P, rng)
            self.rescale.append((np.array(C), np.array(ratio, dtype=np.float64), np.array(res)))
            return res

        def spy_lognorm(C):
            if self.inside:
                self.lognorm.append(np.array(C))
            return real_lognorm(C)

        monkeypatch.setattr(t3d, "build", spy_build)
        monkeypatch.setattr(t3d, "rescale_depth", spy_rescale)
        monkeypatch.setattr(t3d, "lognorm", spy_lognorm)


# ------------------------------------------------------------------------------------------------------ constants


def test_the_line_s_constants_and_strings():
    assert (bk.DEPTH_FORM, bk.DEPTH_BASE, bk.DEPTH_ROW, bk.DEPTH_EXPONENTS, bk.DEPTH_FORMS) == (
        "c2-s2", "s2", "R53", (0.5,), ("global",))
    assert (bk.DEPTH_OPTIONS, bk.DEPTH_FORM_OPTION, bk.DEPTH_FILES) == (OPT, FORM, FILES)
    assert bk.DEPTH_SOURCES == {"allowed-s1": ("allowed", "s1")}
    assert bk.EXIT_DEPTH_CHECK == 7 and len({bk.EXIT_DEPTH_CHECK, bk.EXIT_GUARD_FAILED, bk.EXIT_FACTOR_MISMATCH,
                                             bk.EXIT_REFUSED, bk.EXIT_GATE_FAILED, bk.EXIT_GATE_CRASHED}) == 6
    assert issubclass(bk.DepthCheckFailed, bk.Refused) and not issubclass(bk.DepthCheckFailed, bk.GuardFailed)
    assert bk.C2_S2_SHA256 == "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"
    # the earlier lines are as they were: no form, kind, option string or pin changed
    assert bk.FORMS == ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2",
                        "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2")
    assert bk.KINDS == {"identity": 0.0, "c1": 1.0, "c2": 1.0} and bk.FROZEN_FACTOR == {"source": "mab21l2",
                                                                                         "exponent": 1.0}
    assert (bk.G2_OPTION, bk.C2G_OPTION) == ("`--gata4-exponent 2.0`", "`--form c2-s2 --gata4-exponent 2.0`")
    assert bk.GONLY_OPTIONS == {"F1": "--form c2-s2 --gata4-exponent inf", "F2": "--form c2-s2 --size-scale 0.67",
                                "F3": "--form c2-s2 --gata4-exponent inf --size-scale 0.67"}
    assert not hasattr(bk, "captain_order_names") and not hasattr(bk, "wt_den")


def test_the_pin_is_the_committed_decision_file_and_its_entry_carries_the_step():
    """The pin (here replaced by the synthetic file's) is the committed decision file's sha256, the one row R36 pins;
    its allowed s1 entry carries depth: true and its depth ratios are finite positive numbers (shape only)."""
    src = (ROOT / "experiments" / "t3_decomp" / "build_two_knockout.py").read_text()
    pinned = src.split('DEPTH_DECISION_SHA256 = "', 1)[1].split('"', 1)[0]
    assert pinned == "9440ae67625608112666d5235e6ceaa5dd9f47a4b0bbd30e16943f32f35c5a8c"
    assert pinned == src.split('R36_DECISION_SHA256 = "', 1)[1].split('"', 1)[0]
    assert t3d.sha256(COMMITTED_DECISION) == pinned
    dec = json.loads(COMMITTED_DECISION.read_text())
    assert dec["submissions"]["allowed"]["s1"]["depth"] is True
    pooled = dec["submissions"]["parameters"]["depth_ratio_pooled"]
    assert isinstance(pooled, dict) and pooled
    assert all(isinstance(v, float) and math.isfinite(v) and v > 0 for v in pooled.values())


def test_the_whole_string_match_and_the_containment_hazard():
    """dd's string sits inside ddh's and dg's: containment alone could open dd on either entry; the whole-string
    match does not. No new string is, contains or sits inside an earlier line's."""
    assert OPT["dd"] in OPT["ddh"] and OPT["dd"] in OPT["dg"]
    for a in ARMS:
        cell = f"R53 (...; ruling, no run; `{OPT[a]}` on `{FORM}` only with `--captain-ruling R53`)"
        assert {b for b in ARMS if bk.names_whole(cell, OPT[b])} == {a}
    old = [bk.G2_OPTION.strip("`"), bk.C2G_OPTION.strip("`"), *bk.GONLY_OPTIONS.values(),
           bk._outside_form().R22_OPTION.strip("`"), bk.R36_OPTION.strip("`")]
    for a in OPT.values():
        for b in old:
            assert a != b and f"`{a}`" not in f"`{b}`" and f"`{b}`" not in f"`{a}`", (a, b)


# ------------------------------------------------------------------------------------------------------ the ratios


@pytest.mark.parametrize("arm", ARMS)
def test_the_ratios_by_hand_on_the_reference_s_vocabulary(dw, arm):
    params = {"depth_ratio_pooled": POOLED}
    got, rec = bk.depth_ratios_applied(dw.V, params, **{"ddh": {"exponent": 0.5}, "dg": {"form": "global"}}.get(
        arm, {}))
    want = by_hand(dw.V, POOLED, arm)
    assert np.allclose(got, want, rtol=1e-14, atol=0) and got.shape == (dw.V.L,)
    frozen = bs.depth_ratios(dw.V, params)  # the frozen builder's own reading
    assert rec["frozen"] == [float(x) for x in frozen] and rec["labels"] == list(dw.V.labels)
    assert rec["frozen_sha256"] == hashlib.sha256(np.ascontiguousarray(frozen, dtype=np.float64).tobytes()).hexdigest()
    assert rec["applied_sha256"] == hashlib.sha256(np.ascontiguousarray(got).tobytes()).hexdigest()
    assert rec["labels_sha256"] == hashlib.sha256("\n".join(dw.V.labels).encode()).hexdigest()
    assert rec["n_labels"] == dw.V.L == 11 and rec["n_labels_with_a_pooled_ratio"] == 10
    endo = dw.V.labels.index("Endo")
    assert rec["wt_share_of_labels_left_at_1"] == pytest.approx(float(dw.V.p[endo]), rel=1e-14)
    if arm == "dd":
        assert got.tobytes() == np.asarray(frozen, dtype=np.float64).tobytes()  # the frozen array itself
        assert rec["exponent"] == 1.0 and rec["global_ratio"] is None and rec["form"] == "per label"
    elif arm == "ddh":
        assert got.tobytes() == (np.asarray(frozen, dtype=np.float64) ** 0.5).tobytes()
        assert rec["exponent"] == 0.5 and rec["global_ratio"] is None
    else:
        assert np.all(got == got[0]) and rec["global_ratio"] == got[0] and rec["exponent"] is None
        assert rec["form"] == "global" and rec["rule"] == bk.DEPTH_GLOBAL_RULE


def test_the_global_ratio_is_over_all_labels_with_1_where_none(dw):
    """dg's one ratio is the WT-cell-weighted geometric mean of dd's ratios over all labels, 1.0 where a label has no
    pooled ratio (the cold challenge's fixed rule): over WT E8.75's cells the mean log ratio of dg equals dd's."""
    params = {"depth_ratio_pooled": POOLED}
    dd, _ = bk.depth_ratios_applied(dw.V, params)
    dg, rec = bk.depth_ratios_applied(dw.V, params, form="global")
    n = np.asarray(dw.V.count, dtype=np.float64)
    assert float(n @ np.log(dg)) == pytest.approx(float(n @ np.log(dd)), rel=1e-12)
    # by hand on a small vocabulary: two labels, 3 and 1 cells, ratios 0.25 and 1.0 (left at 1): 0.25 ** 0.75
    v = SimpleNamespace(labels=["a", "b"], L=2, count=np.array([3, 1]), p=np.array([0.75, 0.25]),
                        index=lambda names: np.array([["a", "b"].index(x) for x in names if x in ("a", "b")],
                                                     dtype=int))
    got, r2 = bk.depth_ratios_applied(v, {"depth_ratio_pooled": {"a": 0.25}}, form="global")
    assert np.allclose(got, [0.25 ** 0.75] * 2, rtol=1e-14) and r2["wt_share_of_labels_left_at_1"] == 0.25
    assert not np.allclose(got, 0.25)  # not the mean over the labelled ratios only
    got, _ = bk.depth_ratios_applied(v, {"depth_ratio_pooled": {"a": 0.25, "b": 4.0}}, form="global")
    assert np.allclose(got, [0.25 ** 0.75 * 4.0 ** 0.25] * 2, rtol=1e-14)


def test_the_ratios_refuse_what_the_frozen_reading_refuses_and_bad_values(dw):
    with pytest.raises(bk.Refused, match="not a label of the WT reference"):
        bk.depth_ratios_applied(dw.V, {"depth_ratio_pooled": {**POOLED, "NotALabel": 0.5}})
    for bad in (0.0, -0.5, float("nan"), float("inf")):
        with pytest.raises(bk.Refused, match="finite and positive"):
            bk.depth_ratios_applied(dw.V, {"depth_ratio_pooled": {**POOLED, "vCM1": bad}})
    for e in (1.0, 2.0, 0.25):
        with pytest.raises(bk.Refused, match="not this line's"):
            bk.depth_ratios_applied(dw.V, {"depth_ratio_pooled": POOLED}, exponent=e)
    with pytest.raises(bk.Refused, match="not this line's"):
        bk.depth_ratios_applied(dw.V, {"depth_ratio_pooled": POOLED}, form="local")


# ---------------------------------------------------------------------------------------------- the decision file


def test_the_parameters_are_the_decision_file_s_depth_ratios_only(dw):
    dec = json.loads(dw.dec.read_text())
    params, rec = bk.depth_parameters(dec, dw.dec)
    assert params == {"depth_ratio_pooled": POOLED} and set(params) == {"depth_ratio_pooled"}
    assert rec["decision_sha256"] == rec["pinned_sha256"] == dw.sha and rec["entry"] == "submissions.allowed.s1"
    assert rec["n_labels"] == len(POOLED) and rec["pooled"] == POOLED
    text = "\n".join(f"{k}\t{v!r}" for k, v in POOLED.items())
    assert rec["pooled_sha256"] == hashlib.sha256(text.encode()).hexdigest()


def test_the_decision_file_refusals(dw, tmp_path, monkeypatch):
    dec = json.loads(dw.dec.read_text())
    monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", "0" * 64)
    with pytest.raises(bk.Refused, match="not the pinned one"):
        bk.depth_parameters(dec, dw.dec)
    cases = [("no_depth", {"depth": False}, "does not carry the depth step"),
             ("depth_absent", {"depth": None}, "does not carry the depth step"),
             ("depth_str", {"depth": "true"}, "does not carry the depth step"),
             ("no_params", {"parameters": {"stage_alphas_mean": {}}}, "no submissions.parameters.depth_ratio_pooled"),
             ("empty", {"pooled": {}}, "no submissions.parameters.depth_ratio_pooled"),
             ("zero", {"pooled": {**POOLED, "SE": 0.0}}, "not a finite positive number"),
             ("neg", {"pooled": {**POOLED, "SE": -1.0}}, "not a finite positive number"),
             ("bool", {"pooled": {**POOLED, "SE": True}}, "not a finite positive number"),
             ("text", {"pooled": {**POOLED, "SE": "2.0"}}, "not a finite positive number")]
    for name, over, says in cases:
        p = write_decision(tmp_path / f"{name}.json", pooled=over.get("pooled", POOLED), depth=over.get("depth", True),
                           parameters=over.get("parameters"))
        monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", t3d.sha256(p))
        with pytest.raises(bk.Refused, match=says):
            bk.depth_parameters(json.loads(p.read_text()), p)
    with pytest.raises(bk.Refused, match="only allowed-s1"):
        bk.depth_parameters(dec, dw.dec, "allowed-s2")


# ------------------------------------------------------------------------------------------------- the ruling row


def test_the_real_row_r53_opens_the_three_files_and_no_other_row_does(tmp_path):
    """Row R53's first cell as the repository's ledger holds it, read by code (never printed): the check passes it for
    each of the three files, from the ledger itself and alone in a throwaway ledger; any other ruling row is refused,
    and so is R53's cell with any one of its strings, its opener or its ruling words taken away."""
    first = bk.ledger_first_cell("R53", REPO_LEDGER)
    assert first is not None, "the repository's ledger holds no single row R53"
    for arm in ARMS:
        rec = bk.check_depth_line("c2-s2", captain_ruling="R53", tries=REPO_LEDGER, **KW[arm])
        assert (rec["row"], rec["arm"], rec["file"], rec["option"], rec["form_option"]) == (
            "R53", arm, FILES[arm], f"`{OPT[arm]}`", f"`{FORM}`")
        assert rec["ledger_sha256"] == t3d.sha256(REPO_LEDGER)
    alone = tmp_path / "r53.md"
    alone.write_text(f"# ledger\n| {first}| T3 | x |\n")
    assert bk.ledger_first_cell("R53", alone) == first
    for arm in ARMS:
        assert bk.check_depth_line("c2-s2", captain_ruling="R53", tries=alone, **KW[arm])["arm"] == arm
        for other in ("R22", "R36", "R41", "R52", "R54", "R530", "r53", "R53 ", ""):
            with pytest.raises(bk.Refused, match="only row R53"):
                bk.check_depth_line("c2-s2", captain_ruling=other, tries=alone, **KW[arm])
    cuts = [f"`{FORM}`", "`--captain-ruling R53`", "ruling, no run", "the captain's order", *(f"`{o}`" for o in
                                                                                             OPT.values())]
    for cut in cuts:
        assert cut in first, cut
        p = tmp_path / "cut.md"
        p.write_text(f"# ledger\n| {first.replace(cut, 'x')}| T3 | x |\n")
        opened = set()
        for arm in ARMS:
            try:
                bk.check_depth_line("c2-s2", captain_ruling="R53", tries=p, **KW[arm])
                opened.add(arm)
            except bk.Refused:
                pass
        # without the form, the ruling, the words or the opener no file opens; without one file's string, not that file
        want = {a for a in ARMS if cut != f"`{OPT[a]}`"} if cut.startswith("`--depth") else set()
        assert opened == want, (cut, opened)


@pytest.mark.parametrize("rows,ok", [
    (ROW, set(ARMS)),
    (row("dd"), {"dd"}),
    (row("ddh"), {"ddh"}),  # ddh's entry opens ddh only, not dd, whose string it contains
    (row("dg"), {"dg"}),
    (row("dd", "dg"), {"dd", "dg"}),
    (row("ddh", "dg", note=True), set(ARMS)),  # R53's note names dd's string as one whole string: it opens dd too
    ("", set()),  # no row R53
    (ROW + ROW, set()),  # two rows R53
    (row(*ARMS, opener="the captain's named exception of a day, ruling, no run"), set()),
    (row(*ARMS, opener="the captain's ruling of a day, ruling, no run"), set()),
    (row(*ARMS, opener="the captain's order of a day, a run"), set()),
    (ROW.replace("| R53 (", "| R53b ("), set()),
    (row(*ARMS, form="`--form c2-s1`"), set()),
    (row(*ARMS, form="--form c2-s2"), set()),  # the form without its backticks
    (row(*ARMS, form="`--form c2-s2 --depth-from allowed-s1`"), set()),  # the form only inside another string
    (row(*ARMS, ruling="R53"), set()),
    (row(*ARMS, ruling="`--captain-ruling R54`"), set()),
    (ROW.replace(f"`{OPT['dd']}`", OPT["dd"]), {"ddh", "dg"}),  # dd's without its backticks
    (ROW.replace(f"`{OPT['dd']}`", f"`a`{OPT['dd']}`b`"), {"ddh", "dg"}),  # in the text, not as one whole string
    (ROW.replace(f"`{OPT['ddh']}`", "`--depth-from allowed-s1 --depth-exponent 0.50`"), {"dd", "dg"}),
    (ROW.replace(f"`{OPT['dg']}`", "`--depth-form global`"), {"dd", "ddh"}),
    (row(*ARMS, rid="R58"), set()),  # the entry is row R58's: the command names R53
    ("| R53 (the captain's order of a day, ruling, no run; not a try) | T3 | " + " ".join(f"`{o}`" for o in
                                                                                      OPT.values())
     + f" `{FORM}` `--captain-ruling R53` |\n", set()),  # the strings outside the first cell
])
def test_the_ruling_row_check(tmp_path, rows, ok):
    p = ledger(tmp_path, rows)
    for arm in ARMS:
        if arm in ok:
            rec = bk.check_depth_line("c2-s2", captain_ruling="R53", tries=p, **KW[arm])
            assert rec == {"row": "R53", "arm": arm, "file": FILES[arm], "option": f"`{OPT[arm]}`",
                           "form_option": f"`{FORM}`", "depth_from": "allowed-s1",
                           "exponent": 0.5 if arm == "ddh" else None, "form": "global" if arm == "dg" else "per label",
                           "ledger": "experiments/harness/TRIES.md", "ledger_sha256": t3d.sha256(p)}
        else:
            with pytest.raises(bk.Refused, match="no single row R53"):
                bk.check_depth_line("c2-s2", captain_ruling="R53", tries=p, **KW[arm])


def test_the_row_check_is_build_outside_form_s_plus_the_whole_strings(tmp_path, monkeypatch):
    bo = bk._outside_form()
    seen = []
    real = bo.captain_order_names
    monkeypatch.setattr(bo, "captain_order_names", lambda rid, option, path=None: seen.append((rid, option, path))
                        or real(rid, option, path))
    p = ledger(tmp_path)
    for arm in ARMS:
        seen.clear()
        assert bk.check_depth_line("c2-s2", captain_ruling="R53", tries=p, **KW[arm])["arm"] == arm
        assert seen == [("R53", f"`{OPT[arm]}`", p)]
    odd = ledger(tmp_path, ROW.replace(f"`{OPT['dd']}`", f"`a`{OPT['dd']}`b`"), name="odd.md")
    assert real("R53", f"`{OPT['dd']}`", odd)  # captain_order_names alone passes it; the whole-string match does not
    with pytest.raises(bk.Refused, match="no single row R53"):
        bk.check_depth_line("c2-s2", "allowed-s1", captain_ruling="R53", tries=odd)


# ------------------------------------------------------------------------------------------------------- refusals


def test_without_the_options_nothing_of_the_line_is_read(monkeypatch):
    monkeypatch.setattr(bk, "_outside_form", lambda: pytest.fail("read on a path without the line"))
    for form in [*bk.FORMS, "c1-s1x1.5", "", None]:
        for kw in ({}, {"captain_ruling": "R53"}, {"captain_ruling": "R22"}, {"composition_from": "allowed-s1"},
                   {"gata4_exponent": 2.0}, {"size_scale": 0.67}, {"den_ruling": "R44"}):
            assert bk.check_depth_line(form, **kw) is None and not bk.is_depth_command()


@pytest.mark.parametrize("form,kw,says", [
    ("c2-s2", {"depth_exponent": 0.5, "captain_ruling": "R53"}, "only with --depth-from allowed-s1"),
    ("c2-s2", {"depth_form": "global", "captain_ruling": "R53"}, "only with --depth-from allowed-s1"),
    ("c2-s2", {"depth_from": "allowed-s2", "captain_ruling": "R53"}, "only allowed-s1"),
    ("c2-s2", {"depth_from": "allowed_s1", "captain_ruling": "R53"}, "only allowed-s1"),
    ("c2-s2", {"depth_from": "allowed-s1 ", "captain_ruling": "R53"}, "only allowed-s1"),
    ("c2-s2", {"depth_from": True, "captain_ruling": "R53"}, "only allowed-s1"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": 1.0, "captain_ruling": "R53"}, "only 0.5"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": 2.0, "captain_ruling": "R53"}, "only 0.5"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": 0.4999999, "captain_ruling": "R53"}, "only 0.5"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": float("nan"), "captain_ruling": "R53"}, "only 0.5"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": math.inf, "captain_ruling": "R53"}, "only 0.5"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": "x", "captain_ruling": "R53"}, "only 0.5"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": True, "captain_ruling": "R53"}, "only 0.5"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_form": "Global", "captain_ruling": "R53"}, "only global"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_form": "per-label", "captain_ruling": "R53"}, "only global"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_form": "global", "depth_exponent": 0.5, "captain_ruling": "R53"},
     "not taken together"),
    ("c1-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53"}, "on --form c2-s2 only"),
    ("identity-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53"}, "on --form c2-s2 only"),
    ("c2-s1", {"depth_from": "allowed-s1", "depth_form": "global", "captain_ruling": "R53"}, "on --form c2-s2 only"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53", "composition_from": "allowed-s1"},
     "takes no --composition-from"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53", "gata4_exponent": 2.0}, "takes no --gata4"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53", "gata4_exponent": math.inf}, "takes no --gata4"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53", "size_scale": 0.67}, "takes no --size-scale"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53", "den_ruling": "R44"}, "takes no --den-ruling"),
    ("c2-s2", {"depth_from": "allowed-s1"}, "only with --captain-ruling R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": 0.5}, "only with --captain-ruling R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "r53"}, "only row R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53 "}, "only row R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R22"}, "only row R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R36"}, "only row R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R54"}, "only row R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "captain_ruling": "R53", "rows": row("ddh", "dg")}, "no single row R53"),
    ("c2-s2", {"depth_from": "allowed-s1", "depth_exponent": 0.5, "captain_ruling": "R53", "rows": row("dd")},
     "no single row R53"),
])
def test_the_builder_refuses_and_builds_nothing(world, dw, tmp_path, monkeypatch, form, kw, says):
    for name in ("depth_parameters", "depth_step", "run_c2g_guard", "run_gonly_guards", "run_gata4_figures"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on a refusal"))
    kw = dict(kw)
    t = ledger(tmp_path, kw.pop("rows", ROW))
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match=says):
        build(world, dw.dec, out, form, tries=t, **kw)
    assert nothing_written(out) and hooks_restored()


@pytest.mark.parametrize("over,says", [("pin", "not the pinned one"), ("no_depth", "does not carry the depth step"),
                                       ("no_params", "no submissions.parameters.depth_ratio_pooled"),
                                       ("unknown_label", "not a label of the WT reference")])
def test_the_decision_file_refusals_at_the_build_write_nothing(world, dw, tmp_path, monkeypatch, over, says):
    dec = dw.dec
    if over == "pin":
        monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", "0" * 64)
    else:
        dec = write_decision(tmp_path / f"{over}.json", pooled={**POOLED, "NotALabel": 0.5} if over == "unknown_label"
                             else POOLED, depth=over != "no_depth",
                             parameters={"stage_alphas_mean": {}} if over == "no_params" else None)
        monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", t3d.sha256(dec))
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused, match=says):
        build_arm(world, dw, tmp_path, "dd", out, dec=dec)
    assert nothing_written(out) and hooks_restored()


def test_an_existing_output_is_refused_before_the_decision_file_is_read(world, dw, tmp_path, monkeypatch):
    out = tmp_path / FILES["dd"]
    out.write_bytes(b"x")
    monkeypatch.setattr(bk, "depth_parameters", lambda *a, **k: pytest.fail("read past an existing file"))
    with pytest.raises(bk.Refused, match="exists; pass --overwrite"):
        build_arm(world, dw, tmp_path, "dd", out)
    assert out.read_bytes() == b"x"


def test_the_hooks_refuse_another_slot_another_spec_and_other_ratios(dw):
    """The resolve hook adds the step to the barred s2 only and only to C2's spec; the ratios hook takes the pinned
    decision file's ratios only. Each hook is gone after its block, whatever happened in it."""
    params = {"depth_ratio_pooled": POOLED}
    dec = bs.load_decision(dw.dec, tx.RULES)
    hashes = bs.check_hashes(dec, allow_mismatch=False)
    out: dict = {}
    with bk.depth_step(params, None, None, out):
        spec, info = bs.resolve_spec(dec, tx.RULES, "barred", "s2", hashes)
        assert spec == {**bs.barred_table(dec, tx.RULES)["s2"], "depth": True}
        assert info["depth_source"] == bk.DEPTH_SOURCE and out["n_resolves"] == 1
        for branch, slot in (("barred", "s1"), ("allowed", "s1")):
            with pytest.raises(bk.Refused, match="to the barred s2 only"):
                bs.resolve_spec(dec, tx.RULES, branch, slot, hashes)
        r = bs.depth_ratios(dw.V, {"depth_ratio_pooled": dict(POOLED), "factor_files": {}})
        assert np.array_equal(r, by_hand(dw.V, POOLED, "dd")) and out["n_ratios"] == 1
        for other in ({"depth_ratio_pooled": {**POOLED, "SE": 1.9}}, {"depth_ratio_pooled": {}}, {}, None):
            with pytest.raises(bk.Refused, match="not the pinned decision file's"):
                bs.depth_ratios(dw.V, other)
    assert hooks_restored()
    real = bs.barred_table
    try:
        bs.barred_table = lambda d, rl: {**real(d, rl), "s2": {**real(d, rl)["s2"], "depth": True}}
        with bk.depth_step(params, None, None, {}), pytest.raises(bk.Refused, match="not C2's"):
            bs.resolve_spec(dec, tx.RULES, "barred", "s2", hashes)
    finally:
        bs.barred_table = real
    assert hooks_restored()


# --------------------------------------------------------------------------- every existing path as it was, C2 first


def test_without_the_options_c2_is_its_frozen_path_byte_for_byte(world, dw, tmp_path, monkeypatch):
    """C2's command (no depth option, and each depth option given as None) builds what the frozen path builds: the
    combined loader at c2, x = 1, around bs.build_submission on the barred s2, as the builder did before the line."""
    got: dict = {}
    with bk.combined_loader("c2", 1.0, got):
        frozen = bs.build_submission(dw.dec, "barred", "s2", world.wt, world.board, N, tmp_path / "frozen.h5ad",
                                     run_dir=world.run, seed=SEED, threads=2)
    for name in LINE_FUNCS:
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on C2's command"))
    plain = build(world, dw.dec, tmp_path / "c2.h5ad")
    none = build(world, dw.dec, tmp_path / "c2n.h5ad", depth_from=None, depth_exponent=None, depth_form=None)
    a, b, c = ((tmp_path / f).read_bytes() for f in ("frozen.h5ad", "c2.h5ad", "c2n.h5ad"))
    assert a == b == c and frozen["output"]["sha256"] == plain["output"]["sha256"] == none["output"]["sha256"]
    for s in (plain, none):
        assert s["construction"] == frozen["construction"] and s["spec"] == frozen["spec"]
        assert s["construction"]["depth_ratios"] is None and s["combined_build"]["depth_step"] is None
        assert "depth_line" not in s["combined_build"] and "depth" not in s["spec"]
        assert {k: v for k, v in s["inputs"].items() if k != "factor_note"} == frozen["inputs"]
    assert hooks_restored()


OLD_SETS = [(f, {}) for f in bk.FORMS] + [
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}),
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36"}),
]


@pytest.mark.parametrize("form,kw", OLD_SETS, ids=[f"{f}-{kw.get('captain_ruling', 'plain')}" for f, kw in OLD_SETS])
def test_every_existing_command_never_touches_the_line(world, dw, tmp_path, monkeypatch, form, kw):
    t = ledger(tmp_path)
    extra = dict(kw, tries=t) if kw else {}
    plain = build(world, dw.dec, tmp_path / "plain.h5ad", form, **extra)
    for name in LINE_FUNCS:
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on {form}"))
    trapped = build(world, dw.dec, tmp_path / "trapped.h5ad", form, depth_from=None, depth_exponent=None,
                    depth_form=None, **extra)
    assert (tmp_path / "plain.h5ad").read_bytes() == (tmp_path / "trapped.h5ad").read_bytes()
    assert plain["combined_build"] == trapped["combined_build"] and "depth_line" not in plain["combined_build"]
    assert plain["combined_build"]["depth_step"] is None and hooks_restored()


# ------------------------------------------------------------------------------------------- the files at the build


@pytest.mark.parametrize("arm", ARMS)
def test_each_file_is_c2_with_the_depth_step_and_records_its_line(world, dw, tmp_path, monkeypatch, arm):
    c2 = build(world, dw.dec, tmp_path / "c2.h5ad")
    obs = Observed(monkeypatch)
    out = tmp_path / FILES[arm]
    side = build_arm(world, dw, tmp_path, arm, out)
    want = by_hand(dw.V, POOLED, arm)
    # the frozen builder's own record of the ratios applied (labels whose ratio is not 1.0), as by hand
    assert side["construction"]["depth_ratios"] == pytest.approx(
        {t: float(r) for t, r in zip(dw.V.labels, want, strict=True) if r != 1.0}, rel=1e-14)
    # C2's spec with depth: true added; factor, streams, composition, seed, count and board are C2's
    assert side["spec"] == {**c2["spec"], "depth": True}
    assert side["spec_source"] == {**c2["spec_source"], "depth_source": bk.DEPTH_SOURCE}
    assert side["inputs"]["rng_streams"] == c2["inputs"]["rng_streams"] == ["submission:barred:s2",
                                                                             "submission:barred:s2:assemble"]
    assert side["inputs"]["parameters_source"] == bk.DEPTH_PARAMS_SOURCE == "decision file submissions.parameters"
    assert c2["inputs"]["parameters_source"] == "not needed"
    assert {k: v for k, v in side["inputs"].items() if k != "parameters_source"} == {
        k: v for k, v in c2["inputs"].items() if k != "parameters_source"}
    assert side["construction"]["composition"] == c2["construction"]["composition"]
    assert side["construction"]["factor_exponent"] == 1.0 and side["construction"]["stage_alphas"] is None
    cb, cb2 = side["combined_build"], c2["combined_build"]
    assert {k: v for k, v in cb.items() if k not in ("depth_line", "depth_step")} == {
        k: v for k, v in cb2.items() if k != "depth_step"}
    assert cb["depth_step"] == "combined_build.depth_line" and cb2["depth_step"] is None
    line = cb["depth_line"]
    assert line["file"] == FILES[arm] and line["arm"] == arm and line["line"] == bk.DEPTH_LINE
    assert line["captain_ruling"]["row"] == "R53" and line["captain_ruling"]["option"] == f"`{OPT[arm]}`"
    assert line["exponent"] == {"dd": 1.0, "ddh": 0.5, "dg": None}[arm]
    assert line["form"] == ("global" if arm == "dg" else "per label")
    assert line["spec"] == {"replaced": c2["spec"], "added": {"depth": True}, "source": bk.DEPTH_SOURCE}
    assert np.allclose(line["ratios_applied"], want, rtol=1e-14)
    assert line["ratios_applied_sha256"] == hashlib.sha256(np.asarray(line["ratios_applied"]).tobytes()).hexdigest()
    assert line["ratios_frozen"] == [POOLED.get(t, 1.0) for t in dw.V.labels]
    assert line["decision_file"]["decision_sha256"] == dw.sha and line["decision_file"]["pooled"] == POOLED
    assert line["checks"] == {"pre_depth_counts_are_c2s": True, "rows_are_c2s": True, "labels_are_c2s": True,
                              "sources_are_c2s": True, "ratio_is_the_label_s": True, "one_thinning_draw": True,
                              "passed": True}
    em = line["emptied_by_thinning"]
    assert set(em) == {"stay_unthinned", "empty_before", "empty_after", "handling"}
    assert em["handling"] == bk.DEPTH_EMPTY_RULE and em["empty_before"] == 0 == em["empty_after"]
    assert isinstance(em["stay_unthinned"], int)
    assert line["c2_reference_sha256"] == bk.C2_S2_SHA256
    # the very ratios reached rescale_depth, cell by cell (the file's build only: the C2 twin has no depth step)
    assert len(obs.rescale) == 1
    _, per_cell, _ = obs.rescale[0]
    on_disk = json.loads(out.with_suffix(".json").read_text())["combined_build"]["depth_line"]
    assert on_disk["checks"]["passed"] is True and on_disk["file"] == FILES[arm]
    assert json.loads(Path(f"{out}.stamp.json").read_text())["clean"] is True
    assert set(np.unique(per_cell)) <= set(np.asarray(want)) and hooks_restored()


@pytest.mark.parametrize("arm", ARMS)
def test_each_file_draws_c2_s_cells_changes_only_expression_and_builds_one_sha256(world, dw, tmp_path, arm):
    build(world, dw.dec, tmp_path / "c2.h5ad")
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    sa = build_arm(world, dw, tmp_path, arm, tmp_path / "a" / FILES[arm])
    sb = build_arm(world, dw, tmp_path, arm, tmp_path / "b" / FILES[arm])
    assert (tmp_path / "a" / FILES[arm]).read_bytes() == (tmp_path / "b" / FILES[arm]).read_bytes()
    assert sa["output"]["sha256"] == sb["output"]["sha256"] != t3d.sha256(tmp_path / "c2.h5ad")
    a, c = read(tmp_path / "a" / FILES[arm]), read(tmp_path / "c2.h5ad")
    assert list(a.obs_names) == list(c.obs_names) and a.obs.equals(c.obs) and a.shape == c.shape
    assert list(a.var_names) == list(c.var_names)
    assert np.array_equal(a.obsm["spatial_3D"], c.obsm["spatial_3D"])
    assert not np.array_equal(np.asarray(a.X), np.asarray(c.X))


@pytest.mark.parametrize("arm", ARMS)
def test_the_pre_depth_counts_are_c2_s_final_counts_across_two_builds(world, dw, tmp_path, monkeypatch, arm):
    """C2 built alone (the counts its lognorm receives) and the file built alone (the counts its depth step receives):
    the same array, row by row; the depth step then changes them."""
    obs = Observed(monkeypatch)
    build(world, dw.dec, tmp_path / "c2.h5ad")
    assert len(obs.lognorm) == 1 and not obs.rescale
    c2_counts = obs.lognorm[0]
    obs.lognorm.clear()
    build_arm(world, dw, tmp_path, arm, tmp_path / FILES[arm])
    assert len(obs.rescale) == 1  # the file's depth step, once
    pre, per_cell, post = obs.rescale[0]
    assert pre.dtype == c2_counts.dtype and np.array_equal(pre, c2_counts)
    assert not np.array_equal(post, pre)
    # the builder's twin (C2's cells drawn again inside the build) and the file's own lognorm are both seen here
    assert len(obs.lognorm) == 2 and np.array_equal(obs.lognorm[0], c2_counts)
    assert np.array_equal(obs.lognorm[1], post)


def test_dd_is_the_frozen_builder_s_own_depth_step_on_c2_s_spec(world, dw, tmp_path, monkeypatch):
    """dd byte for byte the frozen builder run on C2's spec with depth: true added by hand (bs.inputs_for then reads
    the decision file's parameters and bs.depth_ratios itself, as the allowed branch's step does)."""
    dd = build_arm(world, dw, tmp_path, "dd", tmp_path / FILES["dd"])
    real = bs.barred_table

    def with_depth(dec, rules):
        tab = real(dec, rules)
        return {**tab, "s2": {**tab["s2"], "depth": True}}

    monkeypatch.setattr(bs, "barred_table", with_depth)
    got: dict = {}
    with bk.combined_loader("c2", 1.0, got):
        frozen = bs.build_submission(dw.dec, "barred", "s2", world.wt, world.board, N, tmp_path / "frozen.h5ad",
                                     run_dir=world.run, seed=SEED, threads=2)
    assert frozen["spec"] == dd["spec"] and frozen["inputs"]["parameters_source"] == bk.DEPTH_PARAMS_SOURCE
    assert (tmp_path / "frozen.h5ad").read_bytes() == (tmp_path / FILES["dd"]).read_bytes()
    assert frozen["construction"]["depth_ratios"] == dd["construction"]["depth_ratios"]


# ------------------------------------------------------------------------------------------ the cells kept unthinned


def test_the_cells_whose_draw_would_be_empty_stay_unthinned_and_are_counted(world, dw, tmp_path, monkeypatch):
    """NCC at ratio 0.004: about half its cells' binomial draws come out empty. The number recorded is the number of
    cells below a ratio of 1 whose final counts are their pre-depth counts (at 0.004 no cell keeps all its counts by
    chance, nor at the other ratios below 1), the frozen thin's rule; no row of the file is empty."""
    monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", dw.tiny_sha)
    obs = Observed(monkeypatch)
    out = tmp_path / FILES["dd"]
    side = build_arm(world, dw, tmp_path, "dd", out, dec=dw.tiny)
    pre, per_cell, post = obs.rescale[0]
    kept = np.all(post == pre, axis=1) & (pre.sum(1) > 0) & (per_cell < 1.0)  # a ratio of 1.0 keeps every count
    assert np.all(per_cell[kept] < 0.01)  # only the tiny label's cells
    em = side["combined_build"]["depth_line"]["emptied_by_thinning"]
    assert em["stay_unthinned"] == int(kept.sum()) > 20 and em["empty_before"] == 0 and em["empty_after"] == 0
    a = read(out)
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    assert np.all(X.sum(1) > 0)  # every row of the file has expression


def test_the_number_is_recorded_on_a_test_board_too(world, dw, tmp_path, monkeypatch):
    """On a board whose counts are hidden the frozen builder still names no cell count of the file (n_cells as a
    label); the line records the number of cells that stay unthinned there as everywhere."""
    from dataclasses import replace

    test_board = replace(world.board, split="test")
    assert test_board.counts_hidden and not world.board.counts_hidden
    monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", dw.tiny_sha)
    obs = Observed(monkeypatch)
    side = build_arm(world, dw, tmp_path, "dd", tmp_path / FILES["dd"], dec=dw.tiny, board=test_board)
    pre, per_cell, post = obs.rescale[0]
    kept = np.all(post == pre, axis=1) & (pre.sum(1) > 0) & (per_cell < 1.0)
    em = side["combined_build"]["depth_line"]["emptied_by_thinning"]
    assert em["stay_unthinned"] == int(kept.sum()) > 0
    assert side["inputs"]["n_cells"] in ("max_cells", "not max_cells")


def test_the_recorder_passes_every_call_to_the_real_generator_unchanged():
    a = np.random.default_rng([7, 1])
    b = np.random.default_rng([7, 1])
    rec = bk._DrawRecorder(a)
    C = np.array([[3, 0, 2], [0, 0, 0], [1, 1, 1], [5, 5, 5]], dtype=np.int64)
    p = np.array([[0.5], [0.5], [0.01], [1.0]])
    got, want = rec.binomial(C, p), b.binomial(C, p)
    assert np.array_equal(got, want) and np.array_equal(rec.binomial_row_sums[0], want.sum(1))
    assert np.array_equal(rec.poisson(np.ones((2, 3))), b.poisson(np.ones((2, 3))))
    assert a.bit_generator.state == b.bit_generator.state
    t = bk._twin_generator(a)
    assert t is not a and np.array_equal(t.integers(0, 1000, 20), a.integers(0, 1000, 20))


# --------------------------------------------------------------------------------------------- the checks' teeth


def _built(rows, lab, source):
    return SimpleNamespace(row=np.asarray(rows), lab=np.asarray(lab), source=np.asarray(source, dtype=object))


def test_the_checks_by_hand():
    pre = np.array([[2, 1], [0, 3], [4, 0]], dtype=np.int32)
    post = np.array([[1, 0], [0, 3], [2, 0]], dtype=np.int32)
    ratios = np.array([0.5, 0.25])
    c2 = _built([5, 6, 7], [0, 1, 0], ["wt"] * 3)
    call = {"pre": pre, "ratio": ratios[[0, 1, 0]], "post": post, "draws": [np.array([1, 0, 2])]}
    ok = bk.depth_checks(c2, pre.copy(), _built([5, 6, 7], [0, 1, 0], ["wt"] * 3), call, ratios)
    assert ok["passed"] and ok["counts"] == {"stay_unthinned": 1, "empty_before": 0, "empty_after": 0}
    for what, args in {
        "pre_depth_counts_are_c2s": (c2, pre + 1, _built([5, 6, 7], [0, 1, 0], ["wt"] * 3), call),
        "rows_are_c2s": (c2, pre, _built([5, 6, 8], [0, 1, 0], ["wt"] * 3), call),
        "labels_are_c2s": (c2, pre, _built([5, 6, 7], [0, 1, 1], ["wt"] * 3),
                           {**call, "ratio": ratios[[0, 1, 1]]}),
        "sources_are_c2s": (c2, pre, _built([5, 6, 7], [0, 1, 0], ["wt", "wt", "stage:x"]), call),
        "ratio_is_the_label_s": (c2, pre, _built([5, 6, 7], [0, 1, 0], ["wt"] * 3),
                                 {**call, "ratio": ratios[[0, 0, 0]]}),
        "one_thinning_draw": (c2, pre, _built([5, 6, 7], [0, 1, 0], ["wt"] * 3), {**call, "draws": []}),
    }.items():
        res = bk.depth_checks(*args, ratios)
        assert res[what] is False and res["passed"] is False, what
    shape = bk.depth_checks(c2, pre[:2], _built([5, 6, 7], [0, 1, 0], ["wt"] * 3), call, ratios)
    assert shape["pre_depth_counts_are_c2s"] is False


def test_a_twin_in_another_state_stops_the_file_with_exit_7_and_writes_nothing(world, dw, tmp_path, monkeypatch):
    monkeypatch.setattr(bk, "_twin_generator", lambda rng: np.random.default_rng(12345))
    out = tmp_path / FILES["ddh"]
    with pytest.raises(bk.DepthCheckFailed, match="pre_depth_counts_are_c2s, rows_are_c2s"):
        build_arm(world, dw, tmp_path, "ddh", out)
    assert nothing_written(out) and hooks_restored()


def test_a_depth_step_that_does_not_receive_c2_s_counts_stops_the_file(world, dw, tmp_path, monkeypatch):
    """A build whose depth step received other counts than C2's final ones (here the file's factor skipped, so the
    step gets the drawn counts) is stopped before anything is written."""
    real, factor = t3d.build, t3d.apply_gene_factor

    def factor_skipped(ctx, n, rng, f, src=None, src_by_label=None, alphas=None, phi=None, w=1.0, ratio=None):
        if ratio is None:
            return real(ctx, n, rng, f, src, src_by_label, alphas, phi, w, None)
        t3d.apply_gene_factor = lambda C, p, P, r: C
        try:
            return real(ctx, n, rng, f, src, src_by_label, alphas, phi, w, ratio)
        finally:
            t3d.apply_gene_factor = factor

    monkeypatch.setattr(t3d, "build", factor_skipped)
    out = tmp_path / FILES["dd"]
    with pytest.raises(bk.DepthCheckFailed, match="pre_depth_counts_are_c2s"):
        build_arm(world, dw, tmp_path, "dd", out)
    assert nothing_written(out) and hooks_restored()


def test_ratios_other_than_the_pinned_file_s_stop_the_build(world, dw, tmp_path, monkeypatch):
    """The ratios the build reads (bs.run_parameters with the run folder) differ from those the pre-check read from
    the pinned decision file (bs.run_parameters without one): refused before any cell is drawn."""
    real = bs.run_parameters
    monkeypatch.setattr(bs, "run_parameters", lambda dec, run_dir: real(dec, run_dir) if run_dir is None else (
        {"depth_ratio_pooled": {**POOLED, "SE": 1.5}}, bk.DEPTH_PARAMS_SOURCE))
    out = tmp_path / FILES["dg"]
    with pytest.raises(bk.Refused, match="not the pinned decision file's"):
        build_arm(world, dw, tmp_path, "dg", out)
    assert nothing_written(out) and hooks_restored()


def test_a_spec_that_already_carries_a_step_is_refused(world, dw, tmp_path, monkeypatch):
    real = bs.barred_table
    monkeypatch.setattr(bs, "barred_table", lambda dec, rules: {**real(dec, rules), "s2": {
        **real(dec, rules)["s2"], "depth": True}})
    out = tmp_path / FILES["dd"]
    with pytest.raises(bk.Refused, match="not C2's"):
        build_arm(world, dw, tmp_path, "dd", out)
    assert nothing_written(out) and hooks_restored()


# --------------------------------------------------------------------------------------------------------- CLI


def test_the_cli(world, dw, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    boards = th.cli_boards(tmp_path)

    def args(*extra, out):
        return ["build", str(dw.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir", str(boards),
                "--run-dir", str(world.run), "--seed", str(SEED), "--threads", "2", *extra, "--out", str(out)]

    runner = CliRunner()
    for extra in (["--form", "c2-s2", "--depth-from", "allowed-s1"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s2", "--captain-ruling", "R53"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R54"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--depth-exponent", "0.5", "--captain-ruling", "R53"],
                  ["--form", "c2-s2", "--depth-form", "global", "--captain-ruling", "R53"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--depth-exponent", "1.0", "--captain-ruling",
                   "R53"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--depth-form", "local", "--captain-ruling", "R53"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--depth-form", "global", "--depth-exponent", "0.5",
                   "--captain-ruling", "R53"],
                  ["--form", "c1-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R53"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R53", "--gata4-exponent",
                   "2.0"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R53", "--composition-from",
                   "allowed-s1"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R53", "--size-scale", "0.67"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R53", "--den-ruling", "R44"],
                  ["--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R53", "--seed", "1000002"]):
        out = tmp_path / "r.h5ad"
        res = runner.invoke(bk.app, args(*extra, out=out))
        assert res.exit_code == 2 and "refused" in res.output, (extra, res.output)
        assert nothing_written(out)
    lines = {"dd": [], "ddh": ["--depth-exponent", "0.5"], "dg": ["--depth-form", "global"]}
    for arm, opts in lines.items():
        ok = tmp_path / FILES[arm]
        res = runner.invoke(bk.app, args("--form", "c2-s2", "--depth-from", "allowed-s1", *opts, "--captain-ruling",
                                         "R53", out=ok))
        assert res.exit_code == 0, res.output
        assert f"the depth line (row R53): {FILES[arm]}" in res.output and "received C2's counts row by row" in (
            res.output)
        line = json.loads(ok.with_suffix(".json").read_text())["combined_build"]["depth_line"]
        assert line["file"] == FILES[arm] and line["captain_ruling"]["arm"] == arm
    # the checks inside the build fail: exit 7, nothing written
    monkeypatch.setattr(bk, "_twin_generator", lambda rng: np.random.default_rng(12345))
    bad = tmp_path / "stop.h5ad"
    res = runner.invoke(bk.app, args("--form", "c2-s2", "--depth-from", "allowed-s1", "--captain-ruling", "R53",
                                     out=bad))
    assert res.exit_code == 7 and "refused (the depth line's checks inside the build" in res.output, res.output
    assert nothing_written(bad) and hooks_restored()
