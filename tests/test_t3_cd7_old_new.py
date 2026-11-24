"""Old against new: experiments/t3_decomp/build_two_knockout.py as it stands against the module at 0470b57 (main before
row R36's option), read with `git show 0470b57:experiments/t3_decomp/build_two_knockout.py`, written to a temporary
folder and loaded beside it, on every existing form and option set, on the synthetic tables and stage of
tests/test_t3_two_knockout_build.py. Each pair gives the same file bytes, the same sidecar but for the fields that name
the run (when and how it ran, its output path, the builder's own sha256), the same stamp, the same refusal messages and
the same command-line exit codes and output. The two routed differences, each refused by both modules here, are named
in their own tests: --captain-ruling R36 (row R36's option) and --gata4-exponent on --form c2-s2 (row R22's third
line); tests/test_t3_c2g_old_new.py compares against the module with row R36's option. Nothing here builds from
released data or reads the repository's ledger.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "0470b57"
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
FORMS = ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2", "c2-s1", "c2-s1x1.5",
         "c2-s1x2.0", "c2-s2")
R22_LEDGER = ("# a throwaway ledger\n| R22 (the captain's order of a day; ruling, no run; the T3 two-knockout builder "
              "is to take `--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`) | T3 | x |\n")
OPTION_SETS = [(f, {}) for f in FORMS] + [("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22"})]
RUN_FIELDS = ("created", "seconds", "argv", "provenance")  # the sidecar's fields that name the run, not the file
STAMP_FIELDS = ("version", "file", "sha256", "commit", "clean", "why", "unclean", "unclean_count", "start_commit")


@pytest.fixture(scope="module")
def old(tmp_path_factory):
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    env |= {"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:{REL}"], capture_output=True, env=env,
                         check=False)
    assert res.returncode == 0, f"{OLD_COMMIT}:{REL} is not in this checkout's history"
    assert b"composition_from" not in res.stdout and b"R36" not in res.stdout  # the module before the option
    d = tmp_path_factory.mktemp("old_module")
    path = d / "build_two_knockout.py"
    path.write_bytes(res.stdout)
    before, frozen = list(sys.path), sys.modules.get("build_submission")
    # the old module imports the frozen builder by name: it gets the very module object the new one uses (another
    # test module may have loaded a second copy under that name since), so both build through one frozen builder
    sys.modules["build_submission"] = bs
    try:
        mod = _load("_build_two_knockout_at_0470b57", path)
    finally:
        sys.path[:] = before  # the old module puts its own folder first on sys.path: taken out again
        if frozen is None:
            sys.modules.pop("build_submission", None)
        else:
            sys.modules["build_submission"] = frozen
    assert mod.bs is bs and mod.SCRIPT_PATH == path.resolve()
    mod.sha256_of_text = hashlib.sha256(res.stdout).hexdigest()
    return mod


@pytest.fixture(autouse=True)
def old_pins(old, world, monkeypatch):
    """The synthetic tables' sha256 in the old module too (the autouse synthetic_pins sets the new module's)."""
    monkeypatch.setattr(old, "TABLE_SHA256", dict(world.pins))


def _ledger(tmp_path: Path) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text(R22_LEDGER)
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


def test_the_old_module_is_main_s_and_its_forms_and_constants_are_the_new_ones(old):
    assert old.FORMS == bk.FORMS == FORMS and old.KINDS == bk.KINDS and old.BASES == bk.BASES
    assert old.REFUSED_FORMS == bk.REFUSED_FORMS and old.FROZEN_SHA256 == bk.FROZEN_SHA256
    assert (old.G2_ROW, old.G2_EXPONENTS, old.G2_FORM, old.G2_OPTION, old.G2_FILE, old.G2_LINE, old.G2_R30) == (
        bk.G2_ROW, bk.G2_EXPONENTS, bk.G2_FORM, bk.G2_OPTION, bk.G2_FILE, bk.G2_LINE, bk.G2_R30)
    assert (old.WHERE, old.LIST_ENTRY, old.FACTOR_NOTE, old.DEN_NOTE, old.SCRIPT) == (
        bk.WHERE, bk.LIST_ENTRY, bk.FACTOR_NOTE, bk.DEN_NOTE, bk.SCRIPT)
    assert old.sha256_of_text != t3d.sha256(ROOT / REL)  # the new module is not the old one


@pytest.mark.parametrize("form,kw", OPTION_SETS, ids=[f"{f}{'-r22' if kw else ''}" for f, kw in OPTION_SETS])
def test_old_and_new_build_the_same_file_and_record(world, old, tmp_path, form, kw):
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    extra = dict(kw, tries=_ledger(tmp_path)) if kw else {}
    so = old.build_form(world.dec, form, world.wt, world.board, N, tmp_path / "old" / "f.h5ad", run_dir=world.run,
                        seed=SEED, threads=2, **extra)
    sn = bk.build_form(world.dec, form, world.wt, world.board, N, tmp_path / "new" / "f.h5ad", run_dir=world.run,
                       seed=SEED, threads=2, **extra)
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
    assert stamp(tmp_path / "new" / "f.h5ad")["clean"] is True
    assert ("gata4_line" in on_new["combined_build"]) is bool(kw) and "composition_line" not in on_new["combined_build"]


REFUSALS = [
    ("c1-s1x1.5", {}), ("c1-s1x2.0", {}), ("c3-s1", {}), ("", {}), (None, {}), (" c1-s1", {}), ("c2-s2 ", {}),
    ("c1-s2", {"gata4_exponent": 3.0, "captain_ruling": "R22"}),
    ("c1-s2", {"gata4_exponent": 2.0}),
    ("c1-s2", {"captain_ruling": "R22"}),
    ("c2-s2", {"captain_ruling": "R22"}),  # the ruling alone on c2-s2: refused as before
    ("c2-s2", {"gata4_exponent": 3.0, "captain_ruling": "R22"}),  # another value on c2-s2: refused as before
    ("c1-s1", {"gata4_exponent": 2.0, "captain_ruling": "R22"}),
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R30"}),
    ("c1-s2", {"gata4_exponent": "x", "captain_ruling": "R22"}),
    ("c1-s2", {"gata4_exponent": 2.0, "captain_ruling": "R22", "ledger": ""}),
    ("c1-s1", {"seed": 0}),
    ("c2-s1", {"seed": 999}),
]


@pytest.mark.parametrize("form,kw", REFUSALS, ids=[f"{f!r}-{i}" for i, (f, _) in enumerate(REFUSALS)])
def test_old_and_new_refuse_alike(world, old, tmp_path, form, kw):
    kw = dict(kw)
    seed = kw.pop("seed", SEED)
    t = tmp_path / "TRIES.md"
    t.write_text(kw.pop("ledger", R22_LEDGER))
    msgs = []
    for mod in (old, bk):
        out = tmp_path / f"{mod.__name__}.h5ad"
        with pytest.raises(bs.Refused) as e:
            mod.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=world.run, seed=seed, threads=2,
                           tries=t, **kw)
        assert not out.exists() and not out.with_suffix(".json").exists()
        msgs.append(str(e.value))
    assert msgs[0] == msgs[1]


def test_old_and_new_refuse_alike_on_a_changed_table_frozen_file_and_existing_output(world, old, tmp_path,
                                                                                      monkeypatch):
    out = {m: tmp_path / f"{m.__name__}.h5ad" for m in (old, bk)}
    for mod in (old, bk):  # a table other than the pinned one
        with monkeypatch.context() as m:
            m.setattr(mod, "TABLE_SHA256", {**world.pins, bk.TABLE_G: "0" * 64})
            with pytest.raises(bs.Refused, match="pinned table"):
                mod.build_form(world.dec, "c2-s2", world.wt, world.board, N, out[mod], run_dir=world.run, seed=SEED,
                               threads=2)
        assert not out[mod].exists()
    changed = tmp_path / "rules.yaml"
    changed.write_bytes(bs.RULES_PATH.read_bytes() + b"\n# changed\n")
    msgs = []
    for mod in (old, bk):  # a changed frozen file
        with monkeypatch.context() as m:
            m.setattr(bs, "RULES_PATH", changed)
            with pytest.raises(bs.Refused) as e:
                mod.build_form(world.dec, "c1-s1", world.wt, world.board, N, out[mod], run_dir=world.run, seed=SEED,
                               threads=2)
            msgs.append(str(e.value))
    assert msgs[0] == msgs[1] and "frozen file" in msgs[0]
    for mod in (old, bk):  # an existing output
        mod.build_form(world.dec, "c1-s1", world.wt, world.board, N, out[mod], run_dir=world.run, seed=SEED,
                       threads=2)
        with pytest.raises(bs.Refused, match="exists"):
            mod.build_form(world.dec, "c1-s1", world.wt, world.board, N, out[mod], run_dir=world.run, seed=SEED,
                           threads=2)
    assert out[old].read_bytes() == out[bk].read_bytes()


CLI_CASES = [[f"--form={f}"] for f in FORMS] + [
    ["--form=c1-s2", "--gata4-exponent=2.0", "--captain-ruling=R22"],
    ["--form=c1-s1x1.5"], ["--form=c3-s1"],
    ["--form=c1-s1", "--seed=1000002"], ["--form=c1-s1", "--n-cells=max_cells"],
    ["--form=c1-s2", "--gata4-exponent=3.0", "--captain-ruling=R22"],
    ["--form=c1-s2", "--gata4-exponent=2.0"], ["--form=c1-s2", "--captain-ruling=R22"],
    ["--form=c2-s2", "--captain-ruling=R22"], ["--form=c2-s2", "--gata4-exponent=3.0", "--captain-ruling=R22"],
    ["--form=c1-s2", "--gata4-exponent=2.0", "--captain-ruling=R30"],
    ["--form=c1-s1", "--exponent=1.5"], ["--form=c1-s1", "--branch=allowed"], ["--form=c1-s1", "--depth"],
]


@pytest.mark.parametrize("extra", CLI_CASES, ids=[" ".join(c) for c in CLI_CASES])
def test_old_and_new_command_lines_exit_and_print_alike(world, old, tmp_path, monkeypatch, extra):
    from typer.testing import CliRunner

    boards = th.cli_boards(tmp_path)
    t = _ledger(tmp_path)
    got = []
    for mod in (old, bk):
        monkeypatch.setattr(mod, "TRIES_FILE", t)
        d = tmp_path / mod.__name__
        d.mkdir()
        args = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
                str(boards), "--run-dir", str(world.run), "--threads", "2", *extra, "--out", str(d / "f.h5ad")]
        res = CliRunner().invoke(mod.app, args)
        got.append((res.exit_code, res.stdout.replace(str(d), "<out>"), (d / "f.h5ad").exists(),
                    (d / "f.h5ad").read_bytes() if (d / "f.h5ad").exists() else None))
    assert got[0] == got[1]
    code, out, built, _ = got[1]
    if code == 0:  # a build: its lines, the path normalised, are compared above
        assert built and out.splitlines()[-1] == "OK: <out>/f.h5ad" and out.startswith("form ")
    else:
        assert not built and code == 2 and ("refused" in out or out == "")  # ours, or the parser's usage error


def test_the_one_routing_difference_is_the_new_ruling_r36(world, old, tmp_path):
    """The new module routes --captain-ruling R36 to row R36's option; the old one refused that ruling as any other.
    Both refuse it without the option, and build nothing."""
    t = _ledger(tmp_path)
    for mod, says in ((old, "only row R22 orders a file of this builder"),
                      (bk, "taken only with --composition-from allowed-s1")):
        out = tmp_path / f"{mod.__name__}.h5ad"
        with pytest.raises(bs.Refused, match=says):
            mod.build_form(world.dec, "c2-s2", world.wt, world.board, N, out, run_dir=world.run, seed=SEED,
                           threads=2, captain_ruling="R36", tries=t)
        assert not out.exists()


@pytest.mark.parametrize("kw", [{"gata4_exponent": 2.0, "captain_ruling": "R22"}, {"gata4_exponent": 2.0}])
def test_the_other_routing_difference_is_row_r22_s_third_line_on_c2_s2(world, old, tmp_path, kw):
    """The new module routes --gata4-exponent on --form c2-s2 to row R22's third line; the old one refused it as a
    form other than c1-s2. Under a ledger whose row R22 names line 2 only, or without the ruling, both refuse it and
    build nothing (the line's own build is compared in tests/test_t3_c2g_build.py)."""
    t = _ledger(tmp_path)
    for mod, says in ((old, "is taken on --form c1-s2 only, not 'c2-s2'"),
                      (bk, "no single row R22" if kw.get("captain_ruling") else "only with --captain-ruling R22")):
        out = tmp_path / f"{mod.__name__}.h5ad"
        with pytest.raises(bs.Refused, match=says):
            mod.build_form(world.dec, "c2-s2", world.wt, world.board, N, out, run_dir=world.run, seed=SEED,
                           threads=2, tries=t, **kw)
        assert not out.exists() and not out.with_suffix(".json").exists()
