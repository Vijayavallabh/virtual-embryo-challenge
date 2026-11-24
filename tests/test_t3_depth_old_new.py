"""Old against new, fourth baseline: experiments/t3_decomp/build_two_knockout.py as it stands against the module at
cbc39c8 (main with row R53, before the depth line's options `--depth-from`, `--depth-exponent` and `--depth-form`),
read with `git show cbc39c8:experiments/t3_decomp/build_two_knockout.py`, written to a temporary folder and loaded
beside it, on the synthetic tables and stage of tests/test_t3_two_knockout_build.py. Every form, line 2's option set
(`--form c1-s2 --gata4-exponent 2.0`), line 3's (`--form c2-s2 --gata4-exponent 2.0`, on a table that passes its guard),
row R36's (cD7) and the options on C2's command (F1, F2, F3, on tables that pass their checks) give the same file bytes
(C2's among them), the same sidecar but for the fields that name the run, the same stamp; the gates g1 and g2 the same
output and exit; the refusals the same messages; the command lines the same exit and output. The routed differences,
the three new options, are named in their own tests. Nothing here builds from released data or reads the repository's
ledger.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "cbc39c8d2796"
REL = "experiments/t3_decomp/build_two_knockout.py"
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES", "GIT_CEILING_DIRECTORIES")


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
POOLED = {"vCM1": 0.6, "OFT/RV-CM": 0.75, "aSHF": 0.5, "pSHF": 0.45, "aPHM": 0.46, "pPHM": 0.55, "NCC": 0.56,
          "SE": 2.0, "Foregut": 0.62, "Peri": 1.4}
FORMS = ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2", "c2-s1", "c2-s1x1.5",
         "c2-s1x2.0", "c2-s2")
LEDGER = (
    "# a throwaway ledger\n"
    "| R22 (the captain's order of a day; ruling, no run; the T3 two-knockout builder is to take "
    "`--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`; by T3's third search line it is also "
    "to take `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`; it is also to take "
    "`--form c2-s2 --gata4-exponent inf`, `--form c2-s2 --size-scale 0.67` and "
    "`--form c2-s2 --gata4-exponent inf --size-scale 0.67`, each only with `--captain-ruling R22`) | T3 | x |\n"
    "| R36 (the captain's named exception of a day; ruling, no run, not a try; the T3 two-knockout builder is to take "
    "`--composition-from allowed-s1` with `--form c2-s2` only with `--captain-ruling R36`) | T3 | y |\n"
    "| R53 (the captain's order of a day, ruling, no run: \"Yes\"; the T3 two-knockout builder is to take "
    "`--depth-from allowed-s1`, `--depth-from allowed-s1 --depth-exponent 0.5` and `--depth-from allowed-s1 "
    "--depth-form global` on `--form c2-s2` only with `--captain-ruling R53`; not a try) | T3 | z |\n")


def alternating(hi: float, lo: float, start: int, stop: int) -> np.ndarray:
    v = np.ones(G)
    v[start:stop:2], v[start + 1:stop:2] = hi, lo
    return v


# tables every check of the F-files passes (tests/test_t3_gonly_build.py's MILD) and line 3's guard (the quarter table)
MILD = (alternating(1.02, 1 / 1.02, 0, 10), alternating(1.02, 1 / 1.02, 10, 20))
GONLY = {"F1": {"gata4_exponent": math.inf}, "F2": {"size_scale": 0.67},
         "F3": {"gata4_exponent": math.inf, "size_scale": 0.67}}
OPTION_SETS = [(f, {}, "world") for f in FORMS] + [
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "world"),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "quarter"),  # line 3, on a table its guard passes
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36"}, "world"),
    *[("c2-s2", {**kw, "captain_ruling": "R22"}, "mild") for kw in GONLY.values()],
]
RUN_FIELDS = ("created", "seconds", "argv", "provenance")
STAMP_FIELDS = ("version", "file", "sha256", "commit", "clean", "why", "unclean", "unclean_count", "start_commit")


def write_decision(path: Path) -> Path:
    """A synthetic decision file: D7m (the barred branch), D7 (row R36's m), the allowed s1 entry (depth: true) and
    submissions.parameters with depth ratios (read by the new module's depth line only)."""
    rules = {"D7m": {"outcome": True, "m": tx.M7M, "n_counted": 1}, "D7": {"outcome": True, "m": L2}}
    entry = {"factor": {"source": "gata4", "exponent": 1.0}, "composition": BLOCK, "depth": True, "stage": False,
             "coverage": False}
    return tx._decision(path, rules=rules, submissions={"allowed": {"s1": entry},
                                                        "parameters": {"depth_ratio_pooled": POOLED}})


@pytest.fixture(scope="module")
def old(tmp_path_factory):
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    env |= {"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:{REL}"], capture_output=True, env=env,
                         check=False)
    assert res.returncode == 0, f"{OLD_COMMIT}:{REL} is not in this checkout's history"
    assert b"c2_gonly_factor" in res.stdout and b"depth_from" not in res.stdout  # the F-files' line, not this one
    d = tmp_path_factory.mktemp("old_module_depth")
    path = d / "build_two_knockout.py"
    path.write_bytes(res.stdout)
    before, frozen = list(sys.path), sys.modules.get("build_submission")
    sys.modules["build_submission"] = bs  # one frozen builder for both modules (see tests/test_t3_cd7_old_new.py)
    try:
        mod = _load("_build_two_knockout_at_cbc39c8", path)
    finally:
        sys.path[:] = before
        if frozen is None:
            sys.modules.pop("build_submission", None)
        else:
            sys.modules["build_submission"] = frozen
    assert mod.bs is bs and mod.SCRIPT_PATH == path.resolve()
    mod.sha256_of_text = hashlib.sha256(res.stdout).hexdigest()
    return mod


@pytest.fixture(scope="module")
def cw(world):
    dec = write_decision(world.d / "gata4_decisions_depth_old_new.json")
    runs = {}
    for key, tables in (("quarter", (th.PHI_M, th.PHI_G ** 0.25)), ("mild", MILD)):
        run = world.d / f"run_{key}_depth_old_new"
        th.write_run(run, world.wt, {bk.TABLE_M: tables[0], bk.TABLE_G: tables[1]})
        runs[key] = SimpleNamespace(run=run, pins={k: t3d.sha256(run / "factors" / f"{k}.csv")
                                                   for k in (bk.TABLE_M, bk.TABLE_G)})
    return SimpleNamespace(dec=dec, sha=t3d.sha256(dec), runs=runs)


@pytest.fixture(autouse=True)
def pins(old, world, cw, monkeypatch):
    """The synthetic tables' and decision file's sha256 in both modules (the autouse synthetic_pins sets the new
    module's tables). The old module, read into a temporary folder, would look for row N4's script beside itself:
    it is pointed at the repository's, as the new module's is."""
    monkeypatch.setattr(old, "TABLE_SHA256", dict(world.pins))
    monkeypatch.setattr(old, "N4_SCRIPT", bk.N4_SCRIPT)
    for mod in (old, bk):
        monkeypatch.setattr(mod, "R36_DECISION_SHA256", cw.sha)
    monkeypatch.setattr(bk, "DEPTH_DECISION_SHA256", cw.sha)


def _ledger(tmp_path: Path, text: str = LEDGER) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text(text)
    return p


def _tables(world, cw, old, monkeypatch, which: str) -> Path:
    """The run folder of the named tables, pinned in both modules (the world's are pinned already)."""
    if which == "world":
        return world.run
    for mod in (old, bk):
        monkeypatch.setattr(mod, "TABLE_SHA256", dict(cw.runs[which].pins))
    return cw.runs[which].run


def comparable(side: dict) -> dict:
    s = json.loads(json.dumps(t3d.json_safe(side)))
    for k in RUN_FIELDS:
        s.pop(k)
    s["output"].pop("h5ad")
    s["combined_build"].pop("script_sha256")
    return s


def stamp(path: Path) -> dict:
    s = json.loads(Path(f"{path}.stamp.json").read_text())
    return {k: s.get(k) for k in STAMP_FIELDS}


def test_the_old_module_s_constants_are_the_new_ones(old):
    for name in ("FORMS", "KINDS", "BASES", "REFUSED_FORMS", "FROZEN_SHA256", "TABLE_M", "TABLE_G", "G2_ROW",
                 "G2_EXPONENTS", "G2_FORM", "G2_OPTION", "G2_FILE", "G2_LINE", "G2_R30", "R36_ROW", "R36_FORM",
                 "R36_OPTION", "R36_NAMES", "R36_FILE", "COMPOSITION_FROM", "R36_LOG2_MAGNITUDE", "R36_RATIO",
                 "R36_RATIO_RTOL", "N4_SETS", "G2_FAILED_LINE", "R36_LINE", "C2G_FORM", "C2G_WEIGHTS", "C2G_OPTION",
                 "C2G_FILE", "C2G_MAX_DEN_SHIFT", "C2G_LINE", "C2G_LIST_ENTRY", "C2G_GUARD_RULE", "C2_S2_SHA256",
                 "GONLY_FORM", "GONLY_SCALES", "GONLY_POINTS", "GONLY_OPTIONS", "GONLY_FILES", "GONLY_RANGE_FACTOR",
                 "GONLY_LIFT_WORDS", "GONLY_LINE", "GONLY_LIST_ENTRY", "GONLY_RANGE_RULE", "GONLY_DEN_RULE",
                 "WHERE", "LIST_ENTRY", "FACTOR_NOTE", "DEN_NOTE", "SCRIPT", "EXIT_GUARD_FAILED",
                 "EXIT_FACTOR_MISMATCH", "EXIT_REFUSED", "EXIT_GATE_FAILED", "EXIT_GATE_CRASHED"):
        assert getattr(old, name) == getattr(bk, name), name
    assert old.sha256_of_text != t3d.sha256(ROOT / REL)  # the new module is not the old one
    assert not hasattr(old, "check_depth_line") and hasattr(bk, "check_depth_line")
    assert old.C2_S2_SHA256 == "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"


@pytest.mark.parametrize("form,kw,tables", OPTION_SETS, ids=[f"{f}-{i}" for i, (f, _, _) in enumerate(OPTION_SETS)])
def test_old_and_new_build_the_same_file_and_record(world, old, cw, tmp_path, monkeypatch, form, kw, tables):
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    run = _tables(world, cw, old, monkeypatch, tables)
    extra = dict(kw, tries=_ledger(tmp_path)) if kw else {}
    so = old.build_form(cw.dec, form, world.wt, world.board, N, tmp_path / "old" / "f.h5ad", run_dir=run, seed=SEED,
                        threads=2, **extra)
    sn = bk.build_form(cw.dec, form, world.wt, world.board, N, tmp_path / "new" / "f.h5ad", run_dir=run, seed=SEED,
                       threads=2, **extra)
    assert (tmp_path / "old" / "f.h5ad").read_bytes() == (tmp_path / "new" / "f.h5ad").read_bytes()
    assert so["output"]["sha256"] == sn["output"]["sha256"]
    assert comparable(so) == comparable(sn)
    on_old = json.loads((tmp_path / "old" / "f.json").read_text())
    on_new = json.loads((tmp_path / "new" / "f.json").read_text())
    assert comparable(on_old) == comparable(on_new)
    assert set(on_new) == set(on_old) and set(on_new["combined_build"]) == set(on_old["combined_build"])
    assert on_old["combined_build"]["script_sha256"] == old.sha256_of_text
    assert on_new["combined_build"]["script_sha256"] == t3d.sha256(ROOT / REL)
    assert stamp(tmp_path / "old" / "f.h5ad") == stamp(tmp_path / "new" / "f.h5ad")
    assert stamp(tmp_path / "new" / "f.h5ad")["clean"] is True and "depth_line" not in on_new["combined_build"]
    assert on_new["combined_build"]["depth_step"] is None and "depth" not in on_new["spec"]


REFUSALS = [
    ("c1-s1x1.5", {}), ("c1-s1x2.0", {}), ("c3-s1", {}), ("", {}), (None, {}), (" c2-s2", {}), ("c2-s2 ", {}),
    ("c1-s2", {"gata4_exponent": 3.0, "captain_ruling": "R22"}), ("c1-s2", {"gata4_exponent": 2.0}),
    ("c1-s2", {"captain_ruling": "R22"}), ("c1-s1", {"gata4_exponent": 2.0, "captain_ruling": "R22"}),
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R30"}),
    ("c1-s2", {"gata4_exponent": math.inf, "captain_ruling": "R22"}),
    ("c2-s2", {"gata4_exponent": -math.inf, "captain_ruling": "R22"}),
    ("c2-s2", {"gata4_exponent": float("nan"), "captain_ruling": "R22"}),
    ("c2-s2", {"captain_ruling": "R22"}), ("c2-s2", {"gata4_exponent": 3.0, "captain_ruling": "R22"}),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22", "tables": "world"}),  # line 3's own guard (exit 6)
    ("c2-s2", {"composition_from": "allowed-s2", "captain_ruling": "R36"}),
    ("c2-s2", {"captain_ruling": "R36"}), ("c1-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36"}),
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R22"}),
    ("c2-s2", {"size_scale": 0.5, "captain_ruling": "R22"}), ("c2-s2", {"size_scale": 0.67}),
    ("c2-s2", {"gata4_exponent": math.inf, "captain_ruling": "R41"}),
    ("c2-s2", {"den_ruling": "R44"}), ("c1-s2", {"size_scale": 0.67, "captain_ruling": "R22"}),
    ("c2-s2", {"captain_ruling": "R53"}),  # the depth line's ruling alone: refused as before
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R53"}), ("c1-s2", {"captain_ruling": "R53"}),
    ("c1-s1", {"seed": 0}),
]


@pytest.mark.parametrize("form,kw", REFUSALS, ids=[f"{f!r}-{i}" for i, (f, _) in enumerate(REFUSALS)])
def test_old_and_new_refuse_alike(world, old, cw, tmp_path, monkeypatch, form, kw):
    kw = dict(kw)
    seed = kw.pop("seed", SEED)
    t = _ledger(tmp_path)
    run = _tables(world, cw, old, monkeypatch, kw.pop("tables", "world"))
    msgs = []
    for mod in (old, bk):
        out = tmp_path / f"{mod.__name__}.h5ad"
        with pytest.raises(bs.Refused) as e:
            mod.build_form(cw.dec, form, world.wt, world.board, N, out, run_dir=run, seed=seed, threads=2, tries=t,
                           **kw)
        assert not out.exists() and not out.with_suffix(".json").exists()
        msgs.append((type(e.value).__name__, str(e.value)))
    assert msgs[0] == msgs[1]


CLI_CASES = [[f"--form={f}"] for f in FORMS] + [
    ["--form=c1-s2", "--gata4-exponent=2.0", "--captain-ruling=R22"],
    ["--form=c2-s2", "--composition-from=allowed-s1", "--captain-ruling=R36"],
    ["--form=c1-s1x1.5"], ["--form=c1-s1", "--seed=1000002"], ["--form=c1-s1", "--n-cells=max_cells"],
    ["--form=c1-s2", "--gata4-exponent=3.0", "--captain-ruling=R22"], ["--form=c2-s2", "--captain-ruling=R22"],
    ["--form=c2-s2", "--gata4-exponent=2.0", "--captain-ruling=R22"],  # line 3 on the world's tables: its guard
    ["--form=c2-s2", "--gata4-exponent=inf", "--captain-ruling=R22"],  # F1 on the world's tables: its range guard
    ["--form=c2-s2", "--size-scale=0.5", "--captain-ruling=R22"], ["--form=c2-s2", "--den-ruling=R44"],
    ["--form=c2-s2", "--captain-ruling=R53"], ["--form=c2-s2", "--composition-from=allowed-s2",
                                               "--captain-ruling=R36"],
]


@pytest.mark.parametrize("extra", CLI_CASES, ids=[" ".join(c) for c in CLI_CASES])
def test_old_and_new_command_lines_exit_and_print_alike(world, old, cw, tmp_path, monkeypatch, extra):
    from typer.testing import CliRunner

    boards = th.cli_boards(tmp_path)
    t = _ledger(tmp_path)
    got = []
    for mod in (old, bk):
        monkeypatch.setattr(mod, "TRIES_FILE", t)
        d = tmp_path / mod.__name__
        d.mkdir()
        args = ["build", str(cw.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
                str(boards), "--run-dir", str(world.run), "--threads", "2", *extra, "--out", str(d / "f.h5ad")]
        res = CliRunner().invoke(mod.app, args)
        got.append((res.exit_code, res.stdout.replace(str(d), "<out>"), (d / "f.h5ad").exists(),
                    (d / "f.h5ad").read_bytes() if (d / "f.h5ad").exists() else None))
    assert got[0] == got[1]
    code, out, built, _ = got[1]
    if code == 0:
        assert built and out.splitlines()[-1] == "OK: <out>/f.h5ad"
    else:
        assert not built and code in (2, 6)


@pytest.mark.parametrize("gate", ["g1", "g2"])
def test_old_and_new_gates_print_and_exit_alike(world, old, cw, tmp_path, gate):
    from typer.testing import CliRunner

    args = [gate, str(cw.dec), "--wt", str(world.wt), "--composition-from", "allowed-s1", "--threads", "2"]
    if gate == "g2":
        args += ["--board", "T3:synthetic", "--boards-dir", str(th.cli_boards(tmp_path)), "--run-dir", str(world.run)]
    got = [(r.exit_code, r.stdout) for r in (CliRunner().invoke(m.app, args) for m in (old, bk))]
    assert got[0] == got[1]
    assert len(got[1][1].splitlines()) == (5 if gate == "g1" else 6) and got[1][0] in (0, 1)


# -------------------------------------------------------------------------------------- the routed differences


@pytest.mark.parametrize("kw", [{"depth_from": "allowed-s1"}, {"depth_exponent": 0.5}, {"depth_form": "global"},
                                {"depth_from": "allowed-s1", "depth_exponent": 0.5},
                                {"depth_from": "allowed-s1", "depth_form": "global"}])
def test_the_new_options_are_the_depth_line_s_and_no_old_command_takes_them(world, old, cw, tmp_path, kw):
    """--depth-from, --depth-exponent and --depth-form: the old module had no such parameters (its build_form a
    TypeError, its command line the parser's usage error); the new one refuses them here without --captain-ruling
    R53, and builds nothing either way."""
    from typer.testing import CliRunner

    t = _ledger(tmp_path)
    out = tmp_path / "x.h5ad"
    with pytest.raises(TypeError, match="depth_"):
        old.build_form(cw.dec, "c2-s2", world.wt, world.board, N, out, run_dir=world.run, seed=SEED, threads=2,
                       tries=t, **kw)
    says = ("only with --captain-ruling R53" if "depth_from" in kw else "only with --depth-from allowed-s1")
    with pytest.raises(bs.Refused, match=says):
        bk.build_form(cw.dec, "c2-s2", world.wt, world.board, N, out, run_dir=world.run, seed=SEED, threads=2,
                      tries=t, **kw)
    assert not out.exists() and not out.with_suffix(".json").exists()
    flags = [f"--{k.replace('_', '-')}={v}" for k, v in kw.items()]
    args = ["build", str(cw.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(th.cli_boards(tmp_path)), "--run-dir", str(world.run), "--threads", "2", "--form=c2-s2", *flags,
            "--out", str(out)]
    res_old, res_new = (CliRunner().invoke(m.app, args) for m in (old, bk))
    assert res_old.exit_code == 2 and res_old.stdout == "" and "--depth-" in res_old.output
    assert res_new.exit_code == 2 and "refused" in res_new.stdout and not out.exists()


def test_the_new_module_builds_the_depth_line_s_files_where_the_old_one_cannot(world, old, cw, tmp_path):
    """With --captain-ruling R53 and row R53 in the ledger, the new module builds dd (C2's cells, the depth step added);
    the old one cannot take the option at all. The file differs from C2's in its bytes, as it must."""
    t = _ledger(tmp_path)
    c2 = bk.build_form(cw.dec, "c2-s2", world.wt, world.board, N, tmp_path / "c2.h5ad", run_dir=world.run,
                       seed=SEED, threads=2)
    dd = bk.build_form(cw.dec, "c2-s2", world.wt, world.board, N, tmp_path / "dd.h5ad", run_dir=world.run, seed=SEED,
                       threads=2, captain_ruling="R53", tries=t, depth_from="allowed-s1")
    assert dd["combined_build"]["depth_line"]["file"] == "T3_c2_s2_dd.h5ad"
    assert dd["spec"] == {**c2["spec"], "depth": True} and dd["output"]["sha256"] != c2["output"]["sha256"]
    with pytest.raises(TypeError, match="depth_from"):
        old.build_form(cw.dec, "c2-s2", world.wt, world.board, N, tmp_path / "o.h5ad", run_dir=world.run, seed=SEED,
                       threads=2, captain_ruling="R53", tries=t, depth_from="allowed-s1")
