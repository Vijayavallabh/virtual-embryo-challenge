"""Row R22's second T3 search line in experiments/t3_decomp/build_two_knockout.py, as row R30 amends it: C1 on the
barred s2 with the Gata4 table's weight raised from 1 to 2 (`--gata4-exponent 2.0`, only with `--captain-ruling R22`
on `--form c1-s2`), its figures recorded and gating nothing, its factor integrity check (exit 5), and the launcher's
third argument g2.

Offline, synthetic and fast: the throwaway world of tests/test_t3_two_knockout_build.py (a small synthetic WT file,
synthetic `transfer_mab21l2_all` and `gata4_pooled` tables, a decision file with the frozen hashes) and small ledger
files written under tmp_path. Nothing here reads the repository's own TRIES.md, builds from released data or pins
bytes built on one machine.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
T3 = ROOT / "experiments" / "t3_decomp"
LAUNCHER = T3 / "build_two_knockout.sh"


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
OPTION_PART = ("by T3's second search line the T3 two-knockout builder is also to take `--gata4-exponent 2.0` on "
               "`--form c1-s2`, only with `--captain-ruling R22`")
T1_PART = "the builder is to take `--feng-fraction-mult 2.5` only with `--captain-ruling R22`"
ROW = f"| R22 (the captain's order of a day: a search; ruling, no run, not a try; {T1_PART}; {OPTION_PART}) | T3 | x |\n"


def ledger(tmp_path: Path, rows: str = ROW, name: str = "TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text(f"# a throwaway ledger\n| R21 (another row) | T2 | y |\n{rows}| R30 (a ruling) | T3 | z |\n")
    return p


def build(world, out: Path, form: str = "c1-s2", **kw) -> dict:
    return bk.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=kw.pop("run_dir", world.run),
                         seed=SEED, threads=2, **kw)


def nothing_written(out: Path) -> bool:
    return not out.exists() and not out.with_suffix(".json").exists() and not Path(f"{out}.stamp.json").exists()


# ------------------------------------------------------------------------------------------------------ constants


def test_the_line_s_constants():
    assert (bk.G2_ROW, bk.G2_EXPONENTS, bk.G2_FORM, bk.G2_OPTION, bk.G2_FILE) == (
        "R22", (2.0,), "c1-s2", "`--gata4-exponent 2.0`", "T3_c1_s2_g2.h5ad")
    assert bk.EXIT_FACTOR_MISMATCH == 5 and bk.EXIT_REFUSED == 2 and bk.G2_OLD_LIMIT == 0.05
    assert bk.TRIES_FILE == ROOT / "experiments" / "harness" / "TRIES.md"
    assert issubclass(bk.FactorMismatch, bk.Refused)
    assert "R30" in bk.G2_R30 and "gates nothing" in bk.G2_R30
    # the existing forms and kinds are as they were: the line adds no form
    assert bk.FORMS == ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2",
                        "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2")
    assert bk.KINDS == {"identity": 0.0, "c1": 1.0, "c2": 1.0}


# ------------------------------------------------------------------------------------------------- the ruling row


@pytest.mark.parametrize("rows,ok", [
    (ROW, True),
    (ROW + ROW, False),  # two rows R22
    ("", False),  # no row R22
    (ROW.replace("| R22 (the captain's order", "| R22 (the captain's ruling"), False),
    (ROW.replace("ruling, no run", "a run"), False),
    (ROW.replace(f"; {OPTION_PART}", ""), False),  # the option absent
    (ROW.replace(f"; {OPTION_PART}) | T3 | x |", ") | T3 | `--gata4-exponent 2.0` with `--captain-ruling R22` |"),
     False),  # the option only outside the first cell
    (ROW.replace("`--gata4-exponent 2.0`", "--gata4-exponent 2.0"), False),  # not in backticks
    (ROW.replace(T1_PART, "the T1 values").replace("only with `--captain-ruling R22`)", "only with R22)"),
     False),  # no `--captain-ruling R22` in the first cell
    (ROW.replace("| R22 (", "| R22b ("), False),
])
def test_the_ruling_row_check(tmp_path, rows, ok):
    p = ledger(tmp_path, rows)
    if ok:
        rec = bk.check_gata4("c1-s2", 2.0, "R22", p)
        assert rec == {"row": "R22", "option": "`--gata4-exponent 2.0`", "gata4_exponent": 2.0,
                       "ledger": "experiments/harness/TRIES.md", "ledger_sha256": t3d.sha256(p)}
    else:
        with pytest.raises(bk.Refused, match="no single row R22"):
            bk.check_gata4("c1-s2", 2.0, "R22", p)


def test_the_row_check_is_build_outside_form_s(tmp_path, monkeypatch):
    """The check is build_outside_form.captain_order_names itself, called with the line's option (reused, not copied)."""
    bo = bk._outside_form()
    seen = []
    monkeypatch.setattr(bo, "captain_order_names", lambda row, option, path=None: seen.append((row, option, path))
                        or True)
    p = tmp_path / "absent.md"
    monkeypatch.setattr(bs, "sha256", lambda path: "x")  # the ledger is absent; only the call is checked
    assert bk.check_gata4("c1-s2", 2.0, "R22", p)["row"] == "R22"
    assert seen == [("R22", "`--gata4-exponent 2.0`", p)]
    assert not hasattr(bk, "captain_order_names") and not hasattr(bk, "wt_den")


def test_a_missing_ledger_is_refused(tmp_path):
    with pytest.raises(bk.Refused, match="no single row R22"):
        bk.check_gata4("c1-s2", 2.0, "R22", tmp_path / "absent.md")


# ------------------------------------------------------------------------------------------------------ refusals


def test_neither_flag_reads_nothing(monkeypatch):
    monkeypatch.setattr(bk, "_outside_form", lambda: pytest.fail("read on a path without the option"))
    for form in bk.FORMS:
        assert bk.check_gata4(form) is None


@pytest.mark.parametrize("value", [1.0, 3.0, 2.5, 0.5, 1.5, 4.0, 0.0, -2.0, float("nan"), float("inf"), "x", "",
                                   True, "2,0"])
def test_any_other_value_is_refused(tmp_path, value):
    with pytest.raises(bk.Refused, match="only 2.0"):
        bk.check_gata4("c1-s2", value, "R22", ledger(tmp_path))


@pytest.mark.parametrize("form", [f for f in ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2",
                                              "c1-s1", "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2")])
def test_any_other_form_is_refused(tmp_path, form):
    with pytest.raises(bk.Refused, match="on --form c1-s2 only"):
        bk.check_gata4(form, 2.0, "R22", ledger(tmp_path))


def test_the_option_without_the_ruling_and_the_ruling_without_the_option(tmp_path):
    with pytest.raises(bk.Refused, match="only with --captain-ruling R22"):
        bk.check_gata4("c1-s2", 2.0, None, ledger(tmp_path))
    for form in ("c1-s2", "c2-s2", "identity-s2"):
        with pytest.raises(bk.Refused, match="taken only with --gata4-exponent 2.0"):
            bk.check_gata4(form, None, "R22", ledger(tmp_path))


@pytest.mark.parametrize("ruling", ["R21", "R30", "r22", "R22 ", " R22", "R22b", "", "K23"])
def test_any_other_ruling_is_refused(tmp_path, ruling):
    for value in (2.0, None):
        with pytest.raises(bk.Refused, match="only row R22"):
            bk.check_gata4("c1-s2", value, ruling, ledger(tmp_path))


@pytest.mark.parametrize("kw", [
    {"gata4_exponent": 3.0, "captain_ruling": "R22"},
    {"gata4_exponent": 2.0},
    {"captain_ruling": "R22"},
    {"gata4_exponent": 2.0, "captain_ruling": "R21"},
    {"gata4_exponent": 2.0, "captain_ruling": "R22", "form": "c2-s2"},
    {"gata4_exponent": 2.0, "captain_ruling": "R22", "form": "c1-s1"},
    {"gata4_exponent": 2.0, "captain_ruling": "R22", "rows": ""},
])
def test_the_builder_refuses_and_builds_nothing(world, tmp_path, kw, monkeypatch):
    monkeypatch.setattr(bk, "run_gata4_figures", lambda *a, **k: pytest.fail("figures computed on a refusal"))
    kw = dict(kw)
    form = kw.pop("form", "c1-s2")
    t = ledger(tmp_path, kw.pop("rows", ROW))
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.Refused):
        build(world, out, form, tries=t, **kw)
    assert nothing_written(out)


# ------------------------------------------------------------------------------------------------------ the factor


def test_the_factor_by_hand_nothing_clipped():
    m = np.array([2.0, 0.5, 1.0, 4.0, 4.0, 0.25, 0.25])
    g = np.array([2.0, 0.5, 4.0, 4.0, 0.25, 0.25, 4.0])
    phi, rec = bk.c1_gata4_weighted(m, g, 2.0)
    assert np.array_equal(phi, [8.0, 0.125, 16.0, 64.0, 0.25, 1.0 / 64.0, 4.0])  # beyond [0.25, 4]: not clipped
    assert rec["gata4_weight"] == 2.0 and rec["mab21l2_exponent"] == 1.0 and rec["c"] == 1.0 and rec["kind"] == "c1"
    # at w = 1 it is C1, combine's own array, bit for bit
    assert np.array_equal(bk.c1_gata4_weighted(m, g, 1.0)[0], bk.combine("c1", 1.0, m, g)[0])
    # the synthetic world's tables
    phi, _ = bk.c1_gata4_weighted(PHI_M, PHI_G, 2.0)
    assert np.array_equal(phi, PHI_M * PHI_G ** 2)
    assert phi.max() == 64.0 and phi.min() == 1.0 / 64.0
    with pytest.raises(bk.Refused, match="same genes"):
        bk.c1_gata4_weighted(m, g[:-1], 2.0)
    with pytest.raises(bk.Refused, match="finite and positive"):
        bk.c1_gata4_weighted(m, np.r_[g[:-1], 0.0], 2.0)


def test_the_figures_by_hand():
    m = np.array([2.0, 0.5, 1.0, 1.0])
    g = np.array([2.0, 1.0, 4.0, 1.0])
    den = np.array([1.0, 2.0, 1.0, 4.0])
    f = bk.gata4_figures(m, g, 2.0, den)
    # totals den . phi / sum den, sum den = 8
    assert f["den_total_phi_m"] == pytest.approx((2 + 1 + 1 + 4) / 8)
    assert f["den_total_c1"] == pytest.approx((4 + 1 + 4 + 4) / 8)
    assert f["den_total_factor"] == pytest.approx((8 + 1 + 16 + 4) / 8)
    assert f["den_shift_c1"] == pytest.approx(13 / 8 - 1)
    assert f["den_shift_factor"] == pytest.approx(29 / 8 - 1)  # far beyond 5%: recorded, gates nothing
    assert f["old_limit"] == {"max_abs_den_shift": 0.05, "status": "lifted by row R30: gates nothing"}
    assert f["gates"] is None and "passed" not in f
    assert f["factor_sha256"] == hashlib.sha256(np.array([8.0, 0.5, 16.0, 1.0]).tobytes()).hexdigest()
    with pytest.raises(bk.Refused, match="den"):
        bk.gata4_figures(m, g, 2.0, den[:-1])


# ------------------------------------------------------------------------------------------------- the line's build


class Spy:
    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


def test_the_line_s_file_is_built_with_its_figures_recorded_whatever_they_are(world, tmp_path, monkeypatch):
    spy = Spy(monkeypatch)
    t = ledger(tmp_path)
    out = tmp_path / bk.G2_FILE
    side = build(world, out, gata4_exponent=2.0, captain_ruling="R22", tries=t)
    assert len(spy.seen) == 1 and np.array_equal(spy.seen[0], PHI_M * PHI_G ** 2)
    assert side["output"]["check"]["ok"] and side["construction"]["factor_exponent"] == 1.0
    on_disk = json.loads(out.with_suffix(".json").read_text())
    cb = on_disk["combined_build"]
    assert cb["form"] == "c1-s2" and cb["combination"]["gata4_weight"] == 2.0
    gl = cb["gata4_line"]
    assert gl["captain_ruling"]["row"] == "R22" and gl["captain_ruling"]["ledger_sha256"] == t3d.sha256(t)
    assert gl["gata4_weight"] == 2.0 and gl["mab21l2_exponent"] == 1.0 and gl["file"] == "T3_c1_s2_g2.h5ad"
    assert gl["line"] == bk.G2_LINE and gl["r30"] == bk.G2_R30 and gl["list_entry"] == bk.G2_LIST_ENTRY
    f = gl["figures"]
    sha = hashlib.sha256(np.ascontiguousarray(PHI_M * PHI_G ** 2).tobytes()).hexdigest()
    assert gl["factor_sha256_applied"] == f["factor_sha256"] == sha
    # den as build_outside_form.wt_den makes it, and the totals from it by hand
    genes, den = bk._outside_form().wt_den(world.wt)
    assert np.array_equal(np.array(f["den"]), den) and f["n_genes"] == len(genes) and "genes" not in f
    s = den.sum()
    assert f["den_total_phi_m"] == pytest.approx(den @ PHI_M / s, rel=1e-12)
    assert f["den_total_c1"] == pytest.approx(den @ (PHI_M * PHI_G) / s, rel=1e-12)
    assert f["den_total_factor"] == pytest.approx(den @ (PHI_M * PHI_G ** 2) / s, rel=1e-12)
    assert f["den_shift_factor"] == pytest.approx(f["den_total_factor"] / f["den_total_phi_m"] - 1, rel=1e-12)
    assert f["den_shift_c1"] == pytest.approx(f["den_total_c1"] / f["den_total_phi_m"] - 1, rel=1e-12)
    assert abs(f["den_shift_factor"]) > 0.05  # a large shift, and the file is built
    assert f["old_limit"]["status"] == "lifted by row R30: gates nothing"
    assert f["tables_sha256"] == {bk.TABLE_M: world.pins[bk.TABLE_M], bk.TABLE_G: world.pins[bk.TABLE_G]}
    assert f["wt_sha256"] == t3d.sha256(world.wt)
    assert json.loads(Path(f"{out}.stamp.json").read_text())["clean"] is True
    assert bs.load_factor.__module__ == "build_submission" and bs.load_factor.__name__ == "load_factor"


def test_a_larger_shift_still_does_not_stop_the_build(world, tmp_path, monkeypatch):
    """Whatever the figures say, the file is built: no constant gates it."""
    real = bk.gata4_figures

    def huge(*a, **k):
        return real(*a, **k) | {"den_shift_factor": 10.0}

    monkeypatch.setattr(bk, "gata4_figures", huge)
    out = tmp_path / "big.h5ad"
    side = build(world, out, gata4_exponent=2.0, captain_ruling="R22", tries=ledger(tmp_path))
    assert out.exists() and side["combined_build"]["gata4_line"]["figures"]["den_shift_factor"] == 10.0


def test_the_line_draws_the_same_cells_as_c1_and_differs_in_expression(world, tmp_path):
    import anndata as ad

    build(world, tmp_path / "c1.h5ad")
    build(world, tmp_path / "g2.h5ad", gata4_exponent=2.0, captain_ruling="R22", tries=ledger(tmp_path))
    a, b = ad.read_h5ad(tmp_path / "c1.h5ad"), ad.read_h5ad(tmp_path / "g2.h5ad")
    assert list(a.obs_names) == list(b.obs_names) and a.obs.equals(b.obs)
    assert np.array_equal(a.obsm["spatial_3D"], b.obsm["spatial_3D"])
    assert not np.array_equal(np.asarray(a.X), np.asarray(b.X))


def test_a_factor_other_than_the_figures_exits_5_and_writes_nothing(world, tmp_path, monkeypatch):
    real = bk.run_gata4_figures
    monkeypatch.setattr(bk, "run_gata4_figures", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    out = tmp_path / "x.h5ad"
    with pytest.raises(bk.FactorMismatch, match="not the one the line's figures"):
        build(world, out, gata4_exponent=2.0, captain_ruling="R22", tries=ledger(tmp_path))
    assert nothing_written(out)
    # other genes than the figures': the same exit
    monkeypatch.setattr(bk, "run_gata4_figures", lambda *a, **k: real(*a, **k) | {"genes": ["Z"]})
    with pytest.raises(bk.FactorMismatch, match="genes"):
        build(world, out, gata4_exponent=2.0, captain_ruling="R22", tries=ledger(tmp_path))
    assert nothing_written(out)
    assert bs.load_factor.__module__ == "build_submission"


def test_the_cli(world, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path))
    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(th.cli_boards(tmp_path)), "--run-dir", str(world.run), "--seed", str(SEED), "--threads", "2"]
    runner = CliRunner()
    for extra in (["--form", "c1-s2", "--gata4-exponent", "3.0", "--captain-ruling", "R22"],
                  ["--form", "c1-s2", "--gata4-exponent", "2.0"],
                  ["--form", "c1-s2", "--captain-ruling", "R22"],
                  ["--form", "c2-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R22"],
                  ["--form", "c1-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R30"],
                  ["--form", "c1-s2", "--gata4-exponent", "two", "--captain-ruling", "R22"]):
        out = tmp_path / "r.h5ad"
        res = runner.invoke(bk.app, [*base, *extra, "--out", str(out)])
        assert res.exit_code == 2, (extra, res.output)
        assert nothing_written(out)
    ok = tmp_path / "ok.h5ad"
    res = runner.invoke(bk.app, [*base, "--form", "c1-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R22",
                                 "--out", str(ok)])
    assert res.exit_code == 0, res.output
    assert json.loads(ok.with_suffix(".json").read_text())["combined_build"]["gata4_line"]["gata4_weight"] == 2.0
    real = bk.run_gata4_figures
    monkeypatch.setattr(bk, "run_gata4_figures", lambda *a, **k: real(*a, **k) | {"factor_sha256": "0" * 64})
    bad = tmp_path / "bad.h5ad"
    res = runner.invoke(bk.app, [*base, "--form", "c1-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R22",
                                 "--out", str(bad)])
    assert res.exit_code == 5, res.output
    assert nothing_written(bad)
    # the ruling row gone: refused
    monkeypatch.setattr(bk, "TRIES_FILE", ledger(tmp_path, "", name="old.md"))
    res = runner.invoke(bk.app, [*base, "--form", "c1-s2", "--gata4-exponent", "2.0", "--captain-ruling", "R22",
                                 "--out", str(tmp_path / "n.h5ad")])
    assert res.exit_code == 2 and nothing_written(tmp_path / "n.h5ad")


# ------------------------------------------------------------------------------- every existing path as it was


@pytest.mark.parametrize("form", list(bk.FORMS))
def test_every_existing_form_never_touches_the_line(world, tmp_path, monkeypatch, form):
    """Without the option the line's code is never reached: each form builds the same bytes and sidecar block with
    the line's functions booby-trapped as without, and the sidecar has no line block."""
    plain = build(world, tmp_path / "plain.h5ad", form)
    for name in ("_outside_form", "run_gata4_figures", "gata4_loader", "c1_gata4_weighted", "gata4_figures"):
        monkeypatch.setattr(bk, name, lambda *a, _n=name, **k: pytest.fail(f"{_n} reached on {form}"))
    trapped = build(world, tmp_path / "trapped.h5ad", form, gata4_exponent=None, captain_ruling=None)
    assert (tmp_path / "plain.h5ad").read_bytes() == (tmp_path / "trapped.h5ad").read_bytes()
    assert plain["combined_build"]["combination"] == trapped["combined_build"]["combination"]
    assert "gata4_line" not in plain["combined_build"] and "gata4_line" not in trapped["combined_build"]


def test_combine_returns_what_it_returns_for_each_kind():
    """combine on synthetic tables, by hand (the existing tests' values, again): unchanged by the line."""
    phi, rec = bk.combine("identity", 1.0, PHI_M, PHI_G)
    assert phi is not None and np.array_equal(phi, PHI_M) and rec["gata4_weight"] == 0.0
    phi, rec = bk.combine("c1", 1.0, PHI_M, PHI_G)
    assert np.array_equal(phi, th.C1_HAND) and rec["gata4_weight"] == 1.0 and rec["c"] == 1.0
    lm, lg = np.log(PHI_M), np.log(PHI_G)
    for x in (1.0, 1.5, 2.0):
        phi, rec = bk.combine("c2", x, PHI_M, PHI_G)
        c = np.std(x * lm) / np.std(x * lm + lg)
        assert rec["c"] == pytest.approx(c, rel=1e-13)
        assert np.array_equal(phi, (PHI_M ** (x * c)) * (PHI_G ** c))
    assert set(rec) == {"kind", "x", "gata4_weight", "mab21l2_exponent", "c", "c_rule", "sd_ln_phi_m", "sd_ln_phi_g",
                        "sd_ln_base", "sd_ln_base_plus_g", "corr_ln_m_ln_g", "base_factor", "factor",
                        "den_weighted_mean", "min", "max", "n_up", "n_down"}
    assert rec["den_weighted_mean"] == bk.DEN_NOTE


# ------------------------------------------------------------------------------------------------------ the launcher


def _sh() -> str:
    return LAUNCHER.read_text()


def _code() -> list[str]:
    return [line.split("#")[0] for line in _sh().splitlines()]


def test_the_launcher_parses_and_deletes_nothing():
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], check=False).returncode == 0
    assert not any(re.search(r"(^|[\s;&|(])rm(\s|$)", line) for line in _code())
    assert not any("unlink" in line or "rmdir" in line for line in _code())


def test_the_launcher_s_g2_strings():
    sh = _sh()
    assert 'G2=${3:-}' in sh
    assert 'if [ -n "$G2" ]; then S=${S}_${G2//[^A-Za-z0-9.]/_}; fi' in sh
    assert 'if [ -n "$G2" ]; then O=${O}_g2; fi' in sh
    assert 'if [ "$G2" != g2 ]; then echo "unknown third argument $G2 (g2 only)" >&2; finish 2; fi' in sh
    assert 'if [ "$BASE" != s2 ]; then echo "g2 is for base s2 only" >&2; finish 2; fi' in sh
    assert "--form c1-s2 --gata4-exponent 2.0 \\\n    --captain-ruling R22" in sh
    assert "G2_FILE=T3_c1_s2_g2.h5ad" in sh
    assert "C1_S2=8718330866429dc089eaffd2f96f66936fe54134eba2d50f8126266e6a3e54b8" in sh
    assert "C2_S2=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7" in sh
    assert "[s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2" in sh
    assert "[s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f" in sh
    assert 'XTESTS=(tests/test_t3_gata4_weight.py)' in sh and (ROOT / "tests" / "test_t3_gata4_weight.py").exists()
    assert "if [ $rc = 5 ]; then finish 55; fi" in sh and "if [ $rc = 4 ]; then finish 56; fi" in sh
    order = [sh.index(s) for s in ("gate s1 20", 'gate "$BASE" 30', 'gate_form c1-s2 "$C1_S2" 60',
                                   'gate_form c2-s2 "$C2_S2" 70', 'one_g2 "$f"', "for form in $CANDS")]
    assert order == sorted(order)
    assert sh.index('if [ -e "$O" ]; then echo "refused: $O exists" >&2; finish 96; fi') < sh.index(
        'echo running > "$S/all.exit"')
    # the two C1/C2 references are the build record's
    rec = (ROOT / "experiments" / "builds" / "results" / "t3two_fb01487.SHA256SUMS").read_text().split()
    names = dict(zip(rec[1::2], rec[0::2], strict=True))
    assert names["c1-s2.h5ad"] in sh and names["c2-s2.h5ad"] in sh


def test_build_output_goes_to_dev_null():
    lines = _sh().splitlines()
    for i, line in enumerate(lines):
        if "build_two_knockout.py build" in line:
            block = " ".join(lines[i:i + 3])
            assert "> /dev/null 2>&1" in block


def _run(tmp_path: Path, *args: str) -> tuple[int, Path]:
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    cwd = tmp_path / "cwd"
    cwd.mkdir(exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env |= {"HOME": str(home), "GIT_CEILING_DIRECTORIES": str(tmp_path), "PATH": "/usr/bin:/bin"}
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=cwd, env=env, capture_output=True, timeout=60,
                       check=False)
    return r.returncode, home / "v" / "scratch"


SHA = "0123456789abcdef0123456789abcdef01234567"


@pytest.mark.parametrize("args,folder", [
    ((SHA, "s2", "g3"), "t3two_0123456_s2_g3"),
    ((SHA, "s2", "G2"), "t3two_0123456_s2_G2"),
    ((SHA, "s1", "g2"), "t3two_0123456_s1_g2"),
    ((SHA, "s1x1.5", "g2"), "t3two_0123456_s1x1.5_g2"),
    ((SHA, "s2", "g2", "extra"), "t3two_0123456_s2_g2"),
])
def test_launcher_refusals_leave_a_refused_record_only(tmp_path, args, folder):
    rc, scratch = _run(tmp_path, *args)
    assert rc == 2
    s = scratch / folder
    assert (s / "refused.exit").read_text().strip() == "2"
    assert not (s / "all.exit").exists() and not (s / "head").exists()
    assert sorted(p.name for p in scratch.iterdir()) == [folder]  # never a folder without the third argument
    assert not (tmp_path / "cwd" / "runs").exists()


def test_launcher_without_the_third_argument_refuses_as_before(tmp_path):
    rc, scratch = _run(tmp_path, SHA, "s3")
    assert rc == 2 and (scratch / "t3two_0123456_s3" / "all.exit").read_text().strip() == "2"
    assert not (scratch / "t3two_0123456_s3" / "refused.exit").exists()
    # outside a repository the HEAD check refuses (97), as before, still before anything is built
    rc, scratch = _run(tmp_path, SHA, "s2")
    assert rc == 97 and (scratch / "t3two_0123456_s2" / "all.exit").read_text().strip() == "97"
    rc, scratch = _run(tmp_path, SHA, "s2", "g2")
    assert rc == 97 and (scratch / "t3two_0123456_s2_g2" / "refused.exit").read_text().strip() == "97"
    assert not (scratch / "t3two_0123456_s2_g2" / "all.exit").exists()
    assert not (tmp_path / "cwd" / "runs").exists()
