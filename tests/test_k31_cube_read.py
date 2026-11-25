"""Try K31's read (experiments/t3_decomp/k31_cube_read.py) on hand-made run folders and synthetic arrays: the PASS bar
of F2 and F6 (mean, side, the second largest null, a tie fails), the merged truth (J4), the guard against K18's
committed barred_form (VOID, the base then counted), INCOMPLETE, G0's STOP, the vertices built, the ledger row's
manifest sha256 and the null's permutation seed. No competition file is read and no arm is scored; the J7 dry run
reads committed results only."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t3_decomp"))
import k31_cube_read as kr  # noqa: E402
import sign_cube as sc  # noqa: E402

from vec import harness  # noqa: E402

C, MACH, SHA = "c" * 40, "dgx", "a" * 64
SEEDS = list(range(10))


def test_pass_needs_a_point_eight_seeds_and_the_second_largest_null():
    d = [1.5] * 8 + [0.5, -0.2]
    assert kr.pass_axis(d, [0.1] * 10)["passes"]
    assert not kr.pass_axis([0.9] * 10, [0.0] * 10)["passes"]  # under one point
    assert not kr.pass_axis([2.0] * 7 + [-0.1] * 3, [0.0] * 10)["passes"]  # seven on the mean's side
    nulls = [0.1] * 8 + [3.0, 1.2]
    r = kr.pass_axis([1.2] * 10, nulls)
    assert r["null_second_largest"] == 1.2 and not r["passes"]  # a tie fails
    assert kr.pass_axis([1.21] * 10, nulls)["passes"]  # above the second largest: exceeds nine of the ten
    neg = kr.pass_axis([-1.5] * 10, [0.2] * 10)
    assert neg["passes"] and neg["mean"] < 0  # the sign never matters for passing
    assert kr.pass_axis([1.5] * 9 + [0.0], [0.0] * 10)["seeds_on_side"] == 9  # a zero seed is on no side


def test_the_vertices_built_follow_the_passing_axes_and_g0():
    g0 = {"axes_left": ["S", "M", "G"], "vertices_left": list(sc.NAMED)}
    assert kr.outcome_vertices(["S"], g0) == (["S", "G"], ["ppn", "npp", "npn"])
    assert kr.outcome_vertices(["M"], {**g0, "axes_left": ["S", "M"]}) == (["M"], ["pnp"])
    assert kr.outcome_vertices(["S", "M"], {**g0, "vertices_left": ["npp", "nnn"]}) == (["S", "M", "G"],
                                                                                        ["npp", "nnn"])


def _ledger(tmp: Path, *lines: str) -> Path:
    p = tmp / "TRIES.md"
    p.write_text("".join(lines))
    return p


def test_row_k31_names_the_script_and_one_manifest_sha256(tmp_path):
    good = f"| K31 (one row, one script {kr.ROW_SCRIPT}; cube tables manifest sha256 `{SHA}`) | T3 | x |\n"
    assert kr.row_k31(_ledger(tmp_path, "| K30 (x) | T3 | y |\n", good)) == SHA
    for bad in (good.replace(kr.ROW_SCRIPT, "`other.py`"), good.replace(SHA, "b" * 63),
                good.replace(")", f"; cube tables manifest sha256 `{'b' * 64}`)", 1)):
        with pytest.raises(kr.ct.Refused):
            kr.row_k31(_ledger(tmp_path, bad))
    with pytest.raises(kr.ct.Refused):
        kr.row_k31(_ledger(tmp_path, good, good))
    with pytest.raises(kr.ct.Refused):
        kr.row_k31(tmp_path / "absent.md")


def _g0(outcome="CONTINUE", axes=("S", "M", "G"), verts=sc.NAMED, commit=C):
    return {"outcome": outcome, "axes_left": list(axes), "vertices_left": list(verts), "manifest_sha256": SHA,
            "line": f"Cube G0: {outcome}: x.", "provenance": {"commit": commit, "tree_clean": True, "machine": MACH}}


def test_read_g0_refuses_another_commit_or_machine(tmp_path):
    p = tmp_path / "g0.json"
    p.write_text(json.dumps(_g0()))
    st = {"commit": C, "machine": MACH}
    assert kr.read_g0(p, st)["outcome"] == "CONTINUE"
    for other in ({"commit": "d" * 40, "machine": MACH}, {"commit": C, "machine": "else"}):
        with pytest.raises(kr.ct.Refused):
            kr.read_g0(p, other)
    with pytest.raises(kr.ct.Refused):
        kr.read_g0(tmp_path / "none.json")


def _refs(tmp: Path, base) -> Path:
    refs = tmp / "refs"
    refs.mkdir()
    for n in ("g1", "g2"):
        rows = [{"seed": s, "scores": {kr.BASE: float(base[s])}} for s in SEEDS]
        (refs / f"plain_{n}.json").write_text(json.dumps({"rows": rows}))
    return refs


def _run_folder(tmp: Path, delta: dict, null: dict, base, axes=("S", "M"), g0=None, commit=C, guard_off=0.0,
                not_made=None) -> Path:
    runs = tmp / "runs"
    runs.mkdir()
    (runs / "g0.json").write_text(json.dumps(g0 or _g0()))
    (runs / "run_started.json").write_text(json.dumps({"commit": C, "tree_clean": True, "machine": MACH,
                                                       "axes": list(axes), "tries": kr.tries_for(axes)}))
    for n in ("g1", "g2"):
        (runs / f"plain_{n}").mkdir()
        for s in SEEDS:
            b = float(base[s]) + guard_off
            scores = {harness.FLOOR: 50.0, kr.BASE: b}
            for j in axes:
                scores[kr.ARM[j]] = b + delta[j][s]
                scores[kr.NULL[j]] = b + null[j][s]
            f = {"seed": s, "task": {"key": kr.TASKS[n].key}, "axes": list(axes), "manifest_sha256": SHA,
                 "rows": [{"seed": s, "scores": scores}], "factors": {"not_made": dict(not_made or {})},
                 "provenance": {"commit": commit if (n, s) == ("g2", 9) else C, "tree_clean": True,
                                "machine": MACH}}
            (runs / f"plain_{n}" / f"seed_{s}.json").write_text(json.dumps(f))
    return runs


@pytest.fixture
def world(tmp_path, monkeypatch):
    base = np.linspace(58.0, 61.0, 10)
    refs = _refs(tmp_path, base)
    monkeypatch.setattr(kr, "K18_SHA256", {f"plain_{n}.json": hashlib.sha256((refs / f"plain_{n}.json")
                                                                               .read_bytes()).hexdigest()
                                           for n in ("g1", "g2")})
    return SimpleNamespace(base=base, refs=refs, tmp=tmp_path)


def test_a_passing_stage_axis_builds_the_cube_over_s_and_g(world):
    d = {"S": [1.6] * 9 + [0.4], "M": [0.3] * 10}
    nl = {"S": [0.2] * 10, "M": [0.4] * 10}
    res = kr.reading(_run_folder(world.tmp, d, nl, world.base), refs=world.refs)
    assert not res["void"] and res["outcome"] == "CONTINUE" and res["passed"] == ["S"]
    assert res["build"] == ["ppn", "npp", "npn"] and res["built_axes"] == ["S", "G"]
    assert res["line"] == "Try K31: CONTINUE (PASS S): the cube is built over S, G: ppn, npp, npn."
    assert res["tries_counted"] == ["cube_np", "null_S", "cube_pn", "null_M"]
    assert res["axes"]["S"]["mean"] == pytest.approx(np.mean(d["S"]))


def test_nothing_passing_reads_stop(world):
    d = {"S": [0.8] * 10, "M": [-0.5] * 10}
    res = kr.reading(_run_folder(world.tmp, d, {"S": [0.0] * 10, "M": [0.0] * 10}, world.base), refs=world.refs)
    assert res["outcome"] == "STOP" and res["line"] == kr.LINE_STOP and res["build"] == []


def test_an_axis_matched_by_its_null_does_not_pass(world):
    d = {"S": [2.0] * 10, "M": [0.0] * 10}
    nl = {"S": [2.5, 2.1] + [0.0] * 8, "M": [0.0] * 10}
    res = kr.reading(_run_folder(world.tmp, d, nl, world.base), refs=world.refs)
    assert res["outcome"] == "STOP" and not res["axes"]["S"]["passes"]


def test_the_guard_voids_the_read_and_counts_the_base(world):
    d = {"S": [2.0] * 10, "M": [2.0] * 10}
    res = kr.reading(_run_folder(world.tmp, d, {"S": [0.0] * 10, "M": [0.0] * 10}, world.base, guard_off=0.2),
                     refs=world.refs)
    assert res["void"] and res["outcome"] == "VOID" and res["build"] == []
    assert res["tries_counted"][-1] == kr.BASE and any(w.startswith("guard:") for w in res["why"])


def test_mixed_commits_void_and_a_missing_seed_is_incomplete(world):
    d = {"S": [2.0] * 10, "M": [2.0] * 10}
    z = {"S": [0.0] * 10, "M": [0.0] * 10}
    res = kr.reading(_run_folder(world.tmp, d, z, world.base, commit="e" * 40), refs=world.refs)
    assert res["void"] and any("different commits" in w for w in res["why"])
    runs = world.tmp / "runs"
    (runs / "plain_g1" / "seed_3.json").rename(runs / "seed_3.aside")
    res = kr.reading(runs, refs=world.refs)
    assert res["incomplete"] and res["missing"] == ["plain_g1/seed_3.json"]


def test_an_arm_not_made_cannot_pass(world):
    d = {"S": [2.0] * 10, "M": [0.0] * 10}
    z = {"S": [0.0] * 10, "M": [0.0] * 10}
    res = kr.reading(_run_folder(world.tmp, d, z, world.base, not_made={"null_S": "k out of range"}),
                     refs=world.refs)
    assert res["outcome"] == "STOP" and not res["axes"]["S"]["made"]


def test_g0_stop_reads_stop_with_no_try(world, tmp_path):
    runs = tmp_path / "stop"
    runs.mkdir()
    (runs / "g0.json").write_text(json.dumps(_g0("STOP", axes=("G",), verts=())))
    res = kr.reading(runs, refs=world.refs)
    assert res["outcome"] == "STOP" and res["tries_counted"] == [] and res["line"].startswith("Try K31: STOP (G0):")


def test_only_the_axes_g0_left_are_read(world):
    d = {"M": [1.5] * 10}
    g0 = _g0(axes=("M", "G"), verts=("ppn", "pnp", "pnn"))
    res = kr.reading(_run_folder(world.tmp, d, {"M": [0.1] * 10}, world.base, axes=("M",), g0=g0), refs=world.refs)
    assert res["outcome"] == "CONTINUE" and res["tries_counted"] == ["cube_pn", "null_M"]
    assert res["build"] == ["ppn", "pnp", "pnn"]


def test_read_never_overwrites_and_keeps_g0s_reading(world, tmp_path):
    d = {"S": [2.0] * 10, "M": [2.0] * 10}
    runs = _run_folder(world.tmp, d, {"S": [0.0] * 10, "M": [0.0] * 10}, world.base)
    (runs / "g0.txt").write_text("Cube G0: CONTINUE: x.\n")
    out = tmp_path / "out"
    a = SimpleNamespace(runs=str(runs), out=str(out))
    assert kr.cmd_read(a) == 0  # a reading (here VOID: the committed K18 files are not this test's pins)
    assert (out / "g0.json").read_bytes() == (runs / "g0.json").read_bytes() and (out / "g0.txt").exists()
    assert kr.cmd_read(a) == kr.EXIT_THERE


def test_the_null_permutation_is_seeded_per_draw_seed():
    rng = np.random.default_rng(3)
    B = rng.uniform(0, 1, (60, 4))
    Q = sc.projector(B)
    q = {"S": Q(rng.normal(size=60)), "M": Q(rng.normal(size=60))}
    basis = rng.uniform(size=60)
    for s in (0, 7):
        z = kr.null_axes("S", q, Q, basis, s)
        want = Q(harness.t3_decile_permutation(q["S"], basis, 20261113 + s))
        assert np.array_equal(z[0], want) and np.array_equal(z[1], q["M"])
    assert not np.array_equal(kr.null_axes("M", q, Q, basis, 0)[1], kr.null_axes("M", q, Q, basis, 1)[1])


def test_the_j7_dry_run_holds_on_k21s_committed_reproduction():
    res = kr.dry_run()
    assert res["passes"] and all(res["k18_pins"].values())
    assert all(v["distance"]["max_abs_diff"] == 0 for v in res["tasks"].values())


ST = {"commit": C, "tree_clean": True, "machine": MACH, "utc": "t"}


def _part_start(out: str, axes) -> None:
    """One part's start, as cmd_run makes it: the record written once and read back (an exception exits 1)."""
    kr.write_started(Path(out), ST, axes)


def test_the_four_parts_start_concurrently_without_a_race(tmp_path):
    """Run k31a's crash: the four parts wrote run_started.json through one shared temporary name, and one part's
    rename removed it from under another. Four processes start together on each of forty fresh folders; none fails,
    one record stands and no temporary name is left."""
    import multiprocessing as mp

    ctx = mp.get_context("fork")
    for trial in range(40):
        out = tmp_path / f"run{trial}"
        out.mkdir()
        procs = [ctx.Process(target=_part_start, args=(str(out), ("S", "M"))) for _ in range(4)]
        for p in procs:
            p.start()
        for p in procs:
            p.join(60)
        assert [p.exitcode for p in procs] == [0, 0, 0, 0], trial
        rec = json.loads((out / "run_started.json").read_text())
        assert rec["commit"] == C and rec["axes"] == ["S", "M"] and rec["tries"] == kr.tries_for(("S", "M"))
        assert not list(out.glob(".run_started.json*"))


def test_one_shared_temporary_name_was_the_race_and_a_per_process_name_is_not(tmp_path, monkeypatch):
    path = tmp_path / "run_started.json"
    tmp = tmp_path / ".run_started.json.part"
    tmp.write_bytes(b"a")
    tmp.write_bytes(b"b")  # two parts wrote the shared name before either renamed it
    tmp.replace(path)
    with pytest.raises(FileNotFoundError):
        tmp.replace(path)  # the second rename has nothing to move: k31a's traceback
    other = tmp_path / "other"
    other.mkdir()
    pid = [101]
    monkeypatch.setattr(kr.os, "getpid", lambda: pid[0])
    assert kr.publish_once(other / "run_started.json", b"a") is True
    pid[0] = 202
    assert kr.publish_once(other / "run_started.json", b"b") is False  # the second writer finds the record
    assert (other / "run_started.json").read_bytes() == b"a" and not list(other.glob(".*.part"))


def test_the_start_record_is_written_once_and_checked(tmp_path):
    rec = kr.write_started(tmp_path, ST, ("S",))
    assert rec["tries"] == ["cube_np", "null_S"] and rec["row"] == "K31"
    assert kr.write_started(tmp_path, {**ST, "utc": "later"}, ("S",))["utc"] == "t"  # read back, not rewritten
    for other in ({**ST, "commit": "d" * 40}, {**ST, "machine": "else"}):
        with pytest.raises(kr.ct.Refused):
            kr.write_started(tmp_path, other, ("S",))
    with pytest.raises(kr.ct.Refused):
        kr.write_started(tmp_path, ST, ("S", "M"))


def test_the_start_command_checks_like_run(tmp_path, monkeypatch):
    monkeypatch.setattr(kr.k18, "stamp", lambda: dict(ST))
    monkeypatch.setattr(kr, "row_k31", lambda path=None: SHA)
    out = tmp_path / "R"
    out.mkdir()
    (out / "g0.json").write_text(json.dumps(_g0()))
    assert kr.cmd_start(SimpleNamespace(out=str(out))) == 0
    assert json.loads((out / "run_started.json").read_text())["axes"] == ["S", "M"]
    assert kr.cmd_start(SimpleNamespace(out=str(out))) == 0  # already this run's
    stop = tmp_path / "stop"
    stop.mkdir()
    (stop / "g0.json").write_text(json.dumps(_g0("STOP")))
    assert kr.cmd_start(SimpleNamespace(out=str(stop))) == kr.EXIT_REFUSED
    assert not (stop / "run_started.json").exists()
    monkeypatch.setattr(kr.k18, "stamp", lambda: {**ST, "tree_clean": False})
    assert kr.cmd_start(SimpleNamespace(out=str(tmp_path / "dirty"))) == kr.EXIT_REFUSED
