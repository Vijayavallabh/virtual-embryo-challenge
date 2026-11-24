"""Old against new, third baseline: experiments/t3_decomp/build_two_knockout.py as it stands against the module at
ed6722a (main with row R22's third T3 line, before the options on C2's command `--gata4-exponent inf` and
`--size-scale 0.67`), read with `git show ed6722a:experiments/t3_decomp/build_two_knockout.py`, written to a temporary
folder and loaded beside it, on the synthetic tables and stage of tests/test_t3_two_knockout_build.py. Every form, line
2's option set (`--form c1-s2 --gata4-exponent 2.0`), line 3's (`--form c2-s2 --gata4-exponent 2.0`, on a table that
passes its guard) and row R36's (cD7) give the same file bytes, the same sidecar but for the fields that name the run,
the same stamp; the gates g1 and g2 the same output and exit; the refusals the same messages; the command lines the
same exit and output. The routed differences, `--gata4-exponent inf` on `--form c2-s2` and any `--size-scale`, are
named in their own tests. Nothing here builds from released data or reads the repository's ledger.
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

import pytest

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "ed6722a"
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
SEED, N = th.SEED, th.N
L2 = 1.3826521560990648
BLOCK = json.loads(json.dumps(bs.composition_spec("central", L2, tx.RULES, "Gata4 replicates (D7)")))
FORMS = ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2", "c2-s1", "c2-s1x1.5",
         "c2-s1x2.0", "c2-s2")
LEDGER = (
    "# a throwaway ledger\n"
    "| R22 (the captain's order of a day; ruling, no run; the T3 two-knockout builder is to take "
    "`--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`; by T3's third search line it is also "
    "to take `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`; it is also to take "
    "`--gata4-exponent inf` and `--form c2-s2 --size-scale 0.67`, each only with `--captain-ruling R22`) | T3 | x |\n"
    "| R36 (the captain's named exception of a day; ruling, no run, not a try; the T3 two-knockout builder is to take "
    "`--composition-from allowed-s1` with `--form c2-s2` only with `--captain-ruling R36`) | T3 | y |\n")
OPTION_SETS = [(f, {}, "world") for f in FORMS] + [
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "world"),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"}, "quarter"),  # line 3, on a table its guard passes
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36"}, "world"),
]
RUN_FIELDS = ("created", "seconds", "argv", "provenance")
STAMP_FIELDS = ("version", "file", "sha256", "commit", "clean", "why", "unclean", "unclean_count", "start_commit")


def write_decision(path: Path, block=BLOCK, d7_m=L2) -> Path:
    """The synthetic decision file of tests/test_t3_cd7_build.py: rule D7 with its m and the allowed s1 entry."""
    rules = {"D7m": {"outcome": True, "m": tx.M7M, "n_counted": 1}, "D7": {"outcome": True, "m": d7_m}}
    entry = {"factor": {"source": "gata4", "exponent": 1.0}, "composition": block, "depth": True, "stage": False,
             "coverage": False}
    return tx._decision(path, rules=rules, submissions={"allowed": {"s1": entry}})


@pytest.fixture(scope="module")
def old(tmp_path_factory):
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    env |= {"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:{REL}"], capture_output=True, env=env,
                         check=False)
    assert res.returncode == 0, f"{OLD_COMMIT}:{REL} is not in this checkout's history"
    assert b"c2_gata4_weighted" in res.stdout and b"size_scale" not in res.stdout  # line 3's, not this line's
    d = tmp_path_factory.mktemp("old_module_gonly")
    path = d / "build_two_knockout.py"
    path.write_bytes(res.stdout)
    before, frozen = list(sys.path), sys.modules.get("build_submission")
    sys.modules["build_submission"] = bs  # one frozen builder for both modules (see tests/test_t3_cd7_old_new.py)
    try:
        mod = _load("_build_two_knockout_at_ed6722a", path)
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
    dec = write_decision(world.d / "gata4_decisions_gonly_old_new.json")
    quarter = world.d / "run_quarter_gonly_old_new"
    th.write_run(quarter, world.wt, {bk.TABLE_M: th.PHI_M, bk.TABLE_G: th.PHI_G ** 0.25})
    qpins = {k: t3d.sha256(quarter / "factors" / f"{k}.csv") for k in (bk.TABLE_M, bk.TABLE_G)}
    return SimpleNamespace(dec=dec, sha=t3d.sha256(dec), quarter=quarter, qpins=qpins)


@pytest.fixture(autouse=True)
def pins(old, world, cw, monkeypatch):
    """The synthetic tables' and decision file's sha256 in both modules (the autouse synthetic_pins sets the new
    module's tables). The old module, read into a temporary folder, would look for row N4's script beside itself:
    it is pointed at the repository's, as the new module's is."""
    monkeypatch.setattr(old, "TABLE_SHA256", dict(world.pins))
    monkeypatch.setattr(old, "N4_SCRIPT", bk.N4_SCRIPT)
    for mod in (old, bk):
        monkeypatch.setattr(mod, "R36_DECISION_SHA256", cw.sha)


def _ledger(tmp_path: Path, text: str = LEDGER) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text(text)
    return p


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
                 "WHERE", "LIST_ENTRY", "FACTOR_NOTE", "DEN_NOTE", "SCRIPT", "EXIT_GUARD_FAILED",
                 "EXIT_FACTOR_MISMATCH", "EXIT_REFUSED"):
        assert getattr(old, name) == getattr(bk, name), name
    assert old.sha256_of_text != t3d.sha256(ROOT / REL)  # the new module is not the old one
    assert not hasattr(old, "check_c2_gonly") and hasattr(bk, "check_c2_gonly")


@pytest.mark.parametrize("form,kw,tables", OPTION_SETS, ids=[f"{f}-{kw.get('captain_ruling', 'plain')}"
                                                             for f, kw, _ in OPTION_SETS])
def test_old_and_new_build_the_same_file_and_record(world, old, cw, tmp_path, monkeypatch, form, kw, tables):
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    run = world.run
    if tables == "quarter":
        run = cw.quarter
        for mod in (old, bk):
            monkeypatch.setattr(mod, "TABLE_SHA256", dict(cw.qpins))
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
    assert stamp(tmp_path / "new" / "f.h5ad")["clean"] is True and "c2_gonly_line" not in on_new["combined_build"]
    assert ("c2_gata4_line" in on_new["combined_build"]) is (tables == "quarter")


REFUSALS = [
    ("c1-s1x1.5", {}), ("c1-s1x2.0", {}), ("c3-s1", {}), ("", {}), (None, {}), (" c2-s2", {}), ("c2-s2 ", {}),
    ("c1-s2", {"gata4_exponent": 3.0, "captain_ruling": "R22"}), ("c1-s2", {"gata4_exponent": 2.0}),
    ("c1-s2", {"captain_ruling": "R22"}), ("c1-s1", {"gata4_exponent": 2.0, "captain_ruling": "R22"}),
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R30"}),
    ("c1-s2", {"gata4_exponent": "x", "captain_ruling": "R22"}),
    ("c1-s2", {"gata4_exponent": math.inf, "captain_ruling": "R22"}),  # inf on line 2's form: refused as before
    ("c2-s1", {"gata4_exponent": math.inf, "captain_ruling": "R22"}),
    ("identity-s2", {"gata4_exponent": math.inf, "captain_ruling": "R22"}),
    ("c2-s2", {"gata4_exponent": -math.inf, "captain_ruling": "R22"}),  # -inf on c2-s2: line 3's refusal, as before
    ("c2-s2", {"gata4_exponent": float("nan"), "captain_ruling": "R22"}),
    ("c2-s2", {"gata4_exponent": 1e308, "captain_ruling": "R22"}),
    ("c2-s2", {"captain_ruling": "R22"}),  # the ruling alone on c2-s2
    ("c2-s2", {"gata4_exponent": 3.0, "captain_ruling": "R22"}),  # another value on c2-s2
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R30"}),  # another ruling on c2-s2
    ("c2-s2", {"gata4_exponent": 2.0}),  # line 3 without the ruling
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22", "ledger": ""}),  # line 3 without its row
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22", "tables": "world"}),  # line 3's own guard (exit 6)
    ("c2-s2", {"composition_from": "allowed-s2", "captain_ruling": "R36"}),
    ("c2-s2", {"captain_ruling": "R36"}), ("c1-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36"}),
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R22"}),
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36", "gata4_exponent": 2.0}),
    ("c2-s2", {"gata4_exponent": 2.0, "captain_ruling": "R36"}),
    ("c1-s1", {"seed": 0}),
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36", "decision": "bad_m"}),
    ("c2-s2", {"composition_from": "allowed-s1", "captain_ruling": "R36", "decision": "bad_states"}),
]


@pytest.mark.parametrize("form,kw", REFUSALS, ids=[f"{f!r}-{i}" for i, (f, _) in enumerate(REFUSALS)])
def test_old_and_new_refuse_alike(world, old, cw, tmp_path, monkeypatch, form, kw):
    kw = dict(kw)
    seed = kw.pop("seed", SEED)
    dec = cw.dec
    which = kw.pop("decision", None)
    if which:
        block = {**BLOCK, "log2_magnitude": 0.7} if which == "bad_m" else {**BLOCK, "depleted_states": BLOCK[
            "depleted_states"][::-1]}
        dec = write_decision(tmp_path / f"{which}.json", block)
        for mod in (old, bk):
            monkeypatch.setattr(mod, "R36_DECISION_SHA256", t3d.sha256(dec))
    t = _ledger(tmp_path, kw.pop("ledger", LEDGER))
    kw.pop("tables", None)  # the world's own tables: line 3's guard fails on them
    msgs = []
    for mod in (old, bk):
        out = tmp_path / f"{mod.__name__}.h5ad"
        with pytest.raises(bs.Refused) as e:
            mod.build_form(dec, form, world.wt, world.board, N, out, run_dir=world.run, seed=seed, threads=2,
                           tries=t, **kw)
        assert not out.exists() and not out.with_suffix(".json").exists()
        msgs.append((type(e.value).__name__, str(e.value)))
    assert msgs[0] == msgs[1]


CLI_CASES = [[f"--form={f}"] for f in FORMS] + [
    ["--form=c1-s2", "--gata4-exponent=2.0", "--captain-ruling=R22"],
    ["--form=c2-s2", "--composition-from=allowed-s1", "--captain-ruling=R36"],
    ["--form=c1-s1x1.5"], ["--form=c1-s1", "--seed=1000002"], ["--form=c1-s1", "--n-cells=max_cells"],
    ["--form=c1-s2", "--gata4-exponent=3.0", "--captain-ruling=R22"], ["--form=c2-s2", "--captain-ruling=R22"],
    ["--form=c2-s2", "--gata4-exponent=3.0", "--captain-ruling=R22"],
    ["--form=c2-s2", "--gata4-exponent=2.0", "--captain-ruling=R30"],
    ["--form=c2-s2", "--gata4-exponent=2.0", "--captain-ruling=R22"],  # line 3 on the world's tables: its guard
    ["--form=c1-s2", "--gata4-exponent=inf", "--captain-ruling=R22"],  # inf on line 2's form: as before
    ["--form=c2-s2", "--gata4-exponent=-inf", "--captain-ruling=R22"],
    ["--form=c2-s2", "--composition-from=allowed-s2", "--captain-ruling=R36"],
    ["--form=c2-s2", "--captain-ruling=R36"], ["--form=c1-s1", "--exponent=1.5"],
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


def test_line_3_s_command_line_builds_alike_on_a_table_its_guard_passes(world, old, cw, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    boards = th.cli_boards(tmp_path)
    t = _ledger(tmp_path)
    got = []
    for mod in (old, bk):
        monkeypatch.setattr(mod, "TRIES_FILE", t)
        monkeypatch.setattr(mod, "TABLE_SHA256", dict(cw.qpins))
        d = tmp_path / mod.__name__
        d.mkdir()
        args = ["build", str(cw.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir", str(boards),
                "--run-dir", str(cw.quarter), "--threads", "2", "--form=c2-s2", "--gata4-exponent=2.0",
                "--captain-ruling=R22", "--out", str(d / "f.h5ad")]
        res = CliRunner().invoke(mod.app, args)
        got.append((res.exit_code, res.stdout.replace(str(d), "<out>"), (d / "f.h5ad").read_bytes()))
    assert got[0] == got[1] and got[1][0] == 0


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


@pytest.mark.parametrize("kw", [{"gata4_exponent": math.inf, "captain_ruling": "R22"}, {"gata4_exponent": math.inf},
                                {"gata4_exponent": "inf", "captain_ruling": "R22"}])
def test_the_one_routed_value_is_inf_on_c2_s2(world, old, cw, tmp_path, kw):
    """--gata4-exponent inf on --form c2-s2: the old module refused it as a value other than 2.0; the new one routes
    it to the line, which refuses it here (a ledger without the line's string, or no ruling) and builds nothing."""
    t = _ledger(tmp_path, LEDGER.replace("`--gata4-exponent inf`", "the weight inf"))
    for mod, says in ((old, "only 2.0"), (bk, "no single row R22" if kw.get("captain_ruling") else
                                          "only with --captain-ruling R22")):
        out = tmp_path / f"{mod.__name__}.h5ad"
        with pytest.raises(bs.Refused, match=says):
            mod.build_form(cw.dec, "c2-s2", world.wt, world.board, N, out, run_dir=world.run, seed=SEED, threads=2,
                           tries=t, **kw)
        assert not out.exists() and not out.with_suffix(".json").exists()


def test_the_new_option_is_the_size_scale(world, old, cw, tmp_path):
    """--size-scale: the old module had no such parameter (its command line no such option, the parser's usage
    error); the new one refuses it here without its ruling, and builds nothing either way."""
    from typer.testing import CliRunner

    t = _ledger(tmp_path)
    out = tmp_path / "x.h5ad"
    with pytest.raises(TypeError, match="size_scale"):
        old.build_form(cw.dec, "c2-s2", world.wt, world.board, N, out, run_dir=world.run, seed=SEED, threads=2,
                       tries=t, size_scale=0.67, captain_ruling="R22")
    with pytest.raises(bs.Refused, match="only with --captain-ruling R22"):
        bk.build_form(cw.dec, "c2-s2", world.wt, world.board, N, out, run_dir=world.run, seed=SEED, threads=2,
                      tries=t, size_scale=0.67)
    assert not out.exists()
    boards = th.cli_boards(tmp_path)
    args = ["build", str(cw.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir", str(boards),
            "--run-dir", str(world.run), "--threads", "2", "--form=c2-s2", "--size-scale=0.67", "--out", str(out)]
    res_old, res_new = (CliRunner().invoke(m.app, args) for m in (old, bk))
    assert res_old.exit_code == 2 and res_old.stdout == "" and "--size-scale" in res_old.output
    assert res_new.exit_code == 2 and "refused" in res_new.stdout and not out.exists()
