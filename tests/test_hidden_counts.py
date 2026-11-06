"""s1 build review (2026-11-01): no tool prints a test board's cell limits or a test file's cell count, and the
board's required keys are checked. A test board's limits can encode the held-out stage's cell count, so they are read
in code only (red-team review 2026-10-31)."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
from conftest import git
from scipy import sparse

from vec.boards import Board
from vec.io import check_submission, make_submission, write_submission
from vec.prep import log1p_cp10k

GENES = [f"G{i}" for i in range(20)]
SECRET_MAX = 4321  # stands for a test board's max_cells: never printed


def _file(tmp_path, n, coords=True, name="pred.h5ad"):
    rng = np.random.default_rng(0)
    X = np.log1p(rng.poisson(2.0, (n, len(GENES))).astype(np.float32))
    a = make_submission(X, GENES, rng.normal(size=(n, 3)).astype(np.float32) if coords else None)
    path = tmp_path / name
    write_submission(a, path)
    return path


def _board(key="T2:heart:test_extrap", split="test_extrap", **kw):
    return Board(key=key, task="T2", genes=GENES, needs_coords=True, split=split, setting="heart", min_cells=10,
                 max_cells=SECRET_MAX, **kw)


def test_check_hides_the_limits_of_a_test_board(tmp_path):
    rep = check_submission(_file(tmp_path, 5), _board())
    assert not rep.ok and any("outside the board's limits" in e for e in rep.errors)
    assert not any(str(SECRET_MAX) in e or " 5 cells" in e for e in rep.errors)
    shown = check_submission(_file(tmp_path, 5), _board(key="T2:heart:val_extrap", split="val_extrap"))
    assert any(str(SECRET_MAX) in e for e in shown.errors)  # validation boards keep their numbers


def test_check_enforces_the_board_required_keys(tmp_path):
    path = _file(tmp_path, 50)
    ok = check_submission(path, _board(obsm_required=["spatial_3D"]))
    assert ok.ok, ok.errors
    rep = check_submission(path, _board(obs_required=["sample"], obsm_required=["spatial_3D", "X_other"]))
    assert "obs['sample'] is required by the board spec" in rep.errors
    assert "obsm['X_other'] is required by the board spec" in rep.errors


def test_vec_check_prints_max_cells_not_the_count(tmp_path, monkeypatch):
    from typer.testing import CliRunner

    from vec import cli

    path = _file(tmp_path, SECRET_MAX)
    monkeypatch.setattr(cli, "get_board", lambda key: _board())
    res = CliRunner().invoke(cli.app, ["check", str(path), "--board", "T2:heart:test_extrap"])
    assert res.exit_code == 0, res.output
    assert "max_cells cells" in res.output and str(SECRET_MAX) not in res.output


def test_vec_boards_hides_test_board_limits(tmp_path, monkeypatch):
    from typer.testing import CliRunner

    from vec import cli

    boards = {"T2:heart:test_extrap": _board(), "T2:heart:val_extrap": _board(key="T2:heart:val_extrap",
                                                                              split="val_extrap")}
    monkeypatch.setattr(cli, "load_boards", lambda: boards)
    res = CliRunner().invoke(cli.app, ["boards"])
    assert res.exit_code == 0, res.output
    assert res.output.count(str(SECRET_MAX)) == 1 and "hidden" in res.output


# ---------------------------------------------------------------------------------------------------------------------
# The T3 builder: seed and count (experiments/t3_decomp/build_submission.py)

T3_DIR = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"


@pytest.fixture(scope="module")
def t3b():
    sys.path.insert(0, str(T3_DIR))
    import build_submission

    return build_submission


def test_t3_count_is_read_from_the_board(t3b):
    b = Board(key="T3:bcatenin", task="T3", genes=GENES, needs_coords=True, split="test", max_cells=SECRET_MAX)
    assert t3b.resolve_n_cells("max", b) == SECRET_MAX
    assert t3b.count_label(SECRET_MAX, b) == "max_cells" and t3b.count_label(7, b) == "not max_cells"
    with pytest.raises(t3b.Refused, match="no max_cells"):
        t3b.resolve_n_cells("max", Board(key="T3:bcatenin", task="T3", genes=GENES, needs_coords=True))
    with pytest.raises(t3b.Refused, match="'max' or an integer"):
        t3b.resolve_n_cells("lots", b)
    syn = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, max_cells=3000)
    assert t3b.resolve_n_cells("250", syn) == 250 and t3b.count_label(250, syn) == 250


@pytest.mark.parametrize("seed", [0, 1, 999_999])
def test_t3_build_refuses_a_seed_below_the_offset(t3b, tmp_path, seed):
    b = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, max_cells=3000)
    with pytest.raises(t3b.Refused, match="below 1000000"):
        t3b.build_submission(tmp_path / "d.json", "allowed", "s1", tmp_path / "wt.h5ad", b, 100,
                             tmp_path / "out.h5ad", seed=seed)


def test_t3_cli_defaults_are_the_offset_seed_and_max_cells(t3b):
    import inspect

    params = inspect.signature(t3b.build).parameters
    assert params["seed"].default.default == 1_000_001 and params["n_cells"].default.default == "max"
    assert inspect.signature(t3b.build_submission).parameters["seed"].default == 1_000_001


# ---------------------------------------------------------------------------------------------------------------------
# experiments/builds/build_test.py: dry runs are named, test builds are not, counts are hidden

BUILDS = Path(__file__).resolve().parents[1] / "experiments" / "builds"


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(BUILDS))
    import build_test

    return build_test


@pytest.mark.parametrize("argv, msg", [
    (["t1", "--board", "T1:val", "--out", "runs/x/T1_s2.h5ad"], "must contain 'dry'"),
    (["t1", "--board", "T2:heart:val_extrap", "--out", "runs/x/T1_dry.h5ad"], "takes only"),
    (["t1", "--out", "runs/x/T1_dry.h5ad"], "must not contain 'dry'"),
    (["heart", "--entry", "s1", "--seed", "5", "--out", "runs/x/h.h5ad"], "below 1000000"),
    (["heart", "--entry", "s2_trim", "--seed", "5", "--out", "runs/x/h.h5ad"], "below 1000000"),
])
def test_build_test_refuses_misnamed_dry_runs_and_low_seeds(bt, argv, msg):
    with pytest.raises(SystemExit) as e:
        bt.main(argv)
    assert msg in str(e.value)


def test_t1_build_holds_cells_without_a_feng_group_by_default(bt):
    """R1 (experiments/harness/TRIES.md; the captain's ruling of 2026-11-05): the T1 build's default keeps the cells
    without a Feng group at their E10.5 expression (no momentum; E9.5 is not loaded), in s1 and in s2 (s1's command
    plus --cm-share); --no-hold-unbracketed restores M1 as built before R1."""
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "runs/x/T1_s1.h5ad"]).hold_unbracketed is True
    assert ap.parse_args(["t1", "--cm-share", "0.60", "--out", "runs/x/T1_s2.h5ad"]).hold_unbracketed is True
    assert ap.parse_args(["t1", "--no-hold-unbracketed", "--out", "runs/x/T1.h5ad"]).hold_unbracketed is False

    def no_load(*args, **kw):
        raise AssertionError("under R1 the build loads no earlier stage")

    assert bt._t1_momentum(True, load=no_load) == {}
    calls = []
    kw = bt._t1_momentum(False, load=lambda path, lineage: calls.append((Path(path).name, lineage)) or "E9.5")
    assert kw == {"start": "E9.5", "damping": bt.harness.B0_DAMPING, "start_source": "E9.5_RNA"}
    assert calls == [("E9.5_RNA.h5ad", False)]


# R2's build flag (experiments/harness/TRIES.md R2): --feng-hygiene, off unless given until the build commit after N12.
# The end-to-end runs below build T1 from synthetic stages through build_t1 itself, with every input monkeypatched: no
# released file, board file, Feng statistic or held-flat list is read.

T1_GENES = [f"g{i}" for i in range(60)] + ["held_a", "held_b"]
T1_HELD = ["held_a", "held_b"]  # stands for the pre-registered held-flat list
T1_SECRET = 1187  # stands for T1:test's max_cells: never printed or stored
# released labels (configs/t1_chamber_classes.yaml), as the N2 test uses them
E105 = {"V-CM": 263, "Endocardium-1": 147, "EndoMT Endocardium": 119, "Early VIC": 131, "Epicardium": 103,
        "BEC": 109, "Blood": 173, "Macrophage": 127, "Hepatocyte": 157}
E95 = {"V-CM": 241, "Endocardium": 151, "Proepicardium": 107, "BEC": 113, "Blood": 139, "Hepatocyte": 181}
FENG_GROUPS = ("CM", "Endocardium", "Endothelium", "Epicardium")
UNION = {"file": "experiments/t1_hygiene/results/n12_hygiene_size/hygiene_genes.csv", "sha256": "ab" * 32,
         "genes": ["g5", "g7", "gone_x", "held_a"], "modules": {"D": 2, "N2": 2},  # gone_x: a name E10.5 lacks
         "sources": {"vdb2017_si": "cd" * 32, "bakken2018_s6": "ef" * 32}}  # a synthetic pinned union


def _t1_stage(sizes: dict, seed: int) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    blocks, labels = [], []
    for lab, n in sizes.items():
        mu = rng.gamma(0.4, 1.0, len(T1_GENES))
        C = rng.poisson(1500.0 * mu / mu.sum() * rng.gamma(10.0, 0.1, (n, 1)))
        C[np.arange(n), rng.integers(0, len(T1_GENES), n)] += 1
        blocks.append(C)
        labels += [lab] * n
    C = sparse.csr_matrix(np.vstack(blocks).astype(np.int32))
    a = ad.AnnData(X=log1p_cp10k(C), layers={"counts": C})
    a.var_names, a.obs_names = T1_GENES, [f"c{seed}_{i}" for i in range(C.shape[0])]
    a.obs["celltype"] = labels
    return a


@pytest.fixture
def t1_env(bt, monkeypatch, own_repo):
    """build_t1's inputs, synthetic; `reads` logs, in order, each E10.5 load and each call for Feng's directions (with
    its gene list), and `held_calls` each held-flat list passed to build_submission."""
    rng = np.random.default_rng(7)
    env = SimpleNamespace(ref=_t1_stage(E105, 2), start=_t1_stage(E95, 1), held_calls=[], reads=[], repo=own_repo,
                          board=Board(key="T1:test", task="T1", genes=list(T1_GENES), needs_coords=False,
                                      min_cells=100, max_cells=T1_SECRET),
                          dirs={g: rng.normal(0.0, 1.0, len(T1_GENES)) for g in FENG_GROUPS},
                          build=bt.t1_scope.build_submission)

    def load(path, **kw):
        assert Path(path).name == "E10.5_RNA.h5ad", path
        env.reads.append("E10.5")
        return env.ref.copy()

    def momentum(hold, load=None):
        return {} if hold else {"start": env.start.copy(), "damping": bt.harness.B0_DAMPING,
                                "start_source": "E9.5_RNA"}

    def directions(genes, *a, **k):  # as feng_directions: NaN where a name is not a (synthetic) Feng gene
        env.reads.append(("Feng", list(genes)))
        idx = [T1_GENES.index(x) if x in T1_GENES else -1 for x in genes]
        return ({g: np.array([v[i] if i >= 0 else np.nan for i in idx]) for g, v in env.dirs.items()},
                {"file": "synthetic"})

    def spy(*a, **k):
        env.held_calls.append(list(k["held_flat"]))
        return env.build(*a, **k)

    monkeypatch.setattr(bt, "load_stage", load)
    monkeypatch.setattr(bt, "_t1_momentum", momentum)
    monkeypatch.setattr(bt, "load_boards", lambda: {"T1:test": env.board})
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", directions)
    monkeypatch.setattr(bt.t1_scope, "held_flat_genes", lambda *a, **k: list(T1_HELD))
    monkeypatch.setattr(bt.t1_scope, "build_submission", spy)
    return env


def _reference_file(bt, env, path: Path, hold: bool = True, held: tuple[str, ...] = tuple(T1_HELD)) -> Path:
    """build_t1's file built directly, with the held-flat list `held`: the pre-registered list (the default) is the file
    build_t1 wrote before R2."""
    ref = env.ref.copy()
    shares = bt.harness.t1_chamber_shares(bt.harness._shallow(ref), "E10.5_RNA")
    groups = bt.feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    momentum = {} if hold else {"start": env.start.copy(), "damping": bt.harness.B0_DAMPING,
                                "start_source": "E9.5_RNA"}
    sub, _ = env.build(ref, env.board, shares, bt.harness.count_for(env.board, ref.n_obs),
                       bt.harness.BUILD_SEED_OFFSET + 1, held_flat=list(held), ref_source="E10.5_RNA",
                       directions=env.dirs, direction_fraction=bt.feng_dirs.bracket_fraction(10.5, 12.5),
                       direction_groups=groups, herd=False, **momentum)
    bt.t1_scope.write_checked(sub, path, env.board)
    return path


def _no_read(*a, **k):
    pytest.fail("nothing is read before the refusal")


def test_feng_hygiene_flag_is_none_unless_given_and_follows_r1(bt):
    ap = bt.build_parser()

    def parse(*argv):
        return ap.parse_args(["t1", *argv, "--out", "runs/x/T1_s1.h5ad"]).feng_hygiene

    assert (parse(), parse("--feng-hygiene"), parse("--no-feng-hygiene")) == (None, True, False)
    assert ap.parse_args(["t1", "--cm-share", "0.60", "--out", "runs/x/T1_s2.h5ad"]).feng_hygiene is None
    assert [bt._feng_hygiene_on(f, h) for f, h in ((None, True), (None, False), (False, True), (False, False))] \
        == [True, False, False, False]  # since the build commit of 2026-11-06 (k): on with R1 unless refused
    assert bt._feng_hygiene_on(True, True) is True
    with pytest.raises(SystemExit, match="refused with --no-hold-unbracketed"):
        bt._feng_hygiene_on(True, False)


def test_hygiene_off_reads_nothing_and_holds_the_pre_registered_list_itself(bt, monkeypatch):
    listed = list(T1_HELD)
    monkeypatch.setattr(bt.t1_scope, "held_flat_genes", lambda *a, **k: listed)
    monkeypatch.setattr(bt.hygiene, "load_union", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    assert bt._t1_union(False) is None  # off: no union, no Feng statistic
    held, rec = bt._t1_held(None, ["Blood", *FENG_GROUPS], T1_GENES)  # off: no guard
    assert held is listed and rec is None


@pytest.mark.parametrize("argv, hold", [(["--no-feng-hygiene"], True), (["--no-hold-unbracketed"], False),
                                        (["--no-hold-unbracketed", "--no-feng-hygiene"], False)])
def test_t1_build_without_the_hygiene_writes_todays_file(bt, t1_env, tmp_path, capsys, argv, hold):
    """--no-feng-hygiene keeps the file of 88b057f byte-identical (M1 with R1, before R2), and --no-hold-unbracketed the
    archived s1 of ae6d6b3; the record gains feng_hygiene None and no longer carries the realised composition."""
    today = _reference_file(bt, t1_env, tmp_path / "today.h5ad", hold)
    capsys.readouterr()
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", *argv, "--out", str(out)])
    assert out.read_bytes() == today.read_bytes()
    assert t1_env.held_calls == [T1_HELD]
    assert t1_env.reads == ["E10.5", ("Feng", T1_GENES)]  # as before R2: no names-only call for Feng's groups
    text = out.with_suffix(".build.json").read_text()
    rec = json.loads(text)
    assert rec["feng_hygiene"] is None and rec["hold_unbracketed"] is hold
    assert rec["n_cells"] == "max_cells" and "realised" not in rec
    assert str(T1_SECRET) not in text + capsys.readouterr().out


@pytest.mark.parametrize("argv", [["--feng-hygiene"], []])  # [] : the default since the build commit of (k)
def test_t1_build_with_the_hygiene_holds_the_pinned_union(bt, t1_env, tmp_path, monkeypatch, argv):
    """--feng-hygiene, and the default with R1, build the same file as build_t1 with the union added to the held-flat
    list; R2's checks run before E10.5 is read (the union, then Feng's group names alone), and the record names the
    union, its sources, the commit and the union genes the build's genes lack."""
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: UNION)
    held = ["g5", "g7", "gone_x", "held_a", "held_b"]  # sorted(pre-registered | union)
    want = _reference_file(bt, t1_env, tmp_path / "want.h5ad", held=held)
    today = _reference_file(bt, t1_env, tmp_path / "today.h5ad")
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", *argv, "--out", str(out)])
    assert t1_env.held_calls == [held]
    assert out.read_bytes() == want.read_bytes() != today.read_bytes()
    assert t1_env.reads == [("Feng", [""]), "E10.5", ("Feng", T1_GENES)]
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["feng_hygiene"] == {"file": UNION["file"], "sha256": UNION["sha256"],
                                   "commit": git(t1_env.repo, "rev-parse", "HEAD"), "genes": 4,
                                   "union_genes_absent_from_build": {"count": 1, "genes": ["gone_x"]},
                                   "modules": UNION["modules"], "sources": UNION["sources"]}


def test_t1_build_refuses_the_hygiene_without_r1_before_reading_anything(bt, t1_env, tmp_path, monkeypatch):
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt, "load_boards", _no_read)
    with pytest.raises(SystemExit, match="refused with --no-hold-unbracketed"):
        bt.main(["t1", "--feng-hygiene", "--no-hold-unbracketed", "--out", str(tmp_path / "T1_s1.h5ad")])
    assert not (tmp_path / "T1_s1.h5ad").exists() and t1_env.reads == []


def test_t1_build_refuses_the_hygiene_while_no_union_sha256_is_pinned_before_reading_anything(bt, t1_env, tmp_path,
                                                                                              monkeypatch):
    cfg = bt.hygiene.load_config()
    monkeypatch.setattr(bt.hygiene, "load_config", lambda *a, **k: dict(cfg, union_sha256=None))
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt, "load_boards", _no_read)
    with pytest.raises(SystemExit, match="pins no union_sha256"):
        bt.main(["t1", "--feng-hygiene", "--out", str(tmp_path / "T1_s1.h5ad")])
    assert not (tmp_path / "T1_s1.h5ad").exists() and t1_env.held_calls == [] and t1_env.reads == []


def test_t1_build_refuses_the_hygiene_while_blood_has_a_feng_direction_before_reading_e105(bt, t1_env, tmp_path,
                                                                                           monkeypatch):
    t1_env.dirs["Blood"] = np.zeros(len(T1_GENES))  # its values are never looked at: a names-only guard
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: UNION)
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt, "load_boards", _no_read)
    with pytest.raises(SystemExit, match="Blood is among Feng's direction groups"):
        bt.main(["t1", "--feng-hygiene", "--out", str(tmp_path / "T1_s1.h5ad")])
    assert not (tmp_path / "T1_s1.h5ad").exists() and t1_env.held_calls == []
    assert t1_env.reads == [("Feng", [""])]  # Feng's group names only
    with pytest.raises(SystemExit, match="Blood is among"):  # and again on the groups the build would use
        bt._t1_held(UNION, ["Blood", *FENG_GROUPS], T1_GENES)


def _feng_stats(path: Path, genes: list[str], cells: dict[str, int]) -> None:
    """Synthetic Feng statistics at `path` (cd1_stats.npz) with the MANIFEST.json the filters write beside it: `cells`
    per label in each half of each end of the bracket."""
    rng = np.random.default_rng(3)
    keys, n, cp = [], [], []
    for label, k in cells.items():
        for stage in ("E9.5", "E14.5"):
            prof = rng.gamma(0.6, 20.0, len(genes))
            for half in "AB":
                keys.append(f"{stage}|LV|{label}|{half}")
                n.append(k)
                cp.append(k * prof)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = np.array(n, dtype=np.float64)
    np.savez(path, genes=np.array(genes), keys=np.array(keys), n=n, S_cp=np.array(cp),
             S_log=np.zeros((len(keys), len(genes))), S_nz=np.zeros((len(keys), len(genes))), depth_sum=n * 1000.0)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    (path.parent / "MANIFEST.json").write_text(json.dumps({"outputs": {path.name: {"sha256": sha}}}))


def test_feng_group_names_need_no_t1_file_and_are_the_builds_groups(bt, tmp_path, monkeypatch):
    """The Blood guard's names, from the real feng_dirs.feng_directions on synthetic Feng statistics: its placeholder
    call gives the groups a call on gene names gives (an empty list would fail in dirlib.sym_index)."""
    monkeypatch.setattr(bt.feng_dirs, "REPO", tmp_path)
    genes = [f"f{i}" for i in range(40)]
    cells = {"ventricular_cm": 30, "atrial_cm": 30, "endocardial_ec": 30, "vascular_ec": 30, "epicardial": 30,
             "blood": 2}  # blood below min_cells_ext: no direction
    stats = bt.feng_dirs.clean_stats_path()
    _feng_stats(stats, genes, cells)
    with pytest.raises(IndexError):
        bt.feng_dirs.feng_directions([])
    names = bt._feng_group_names()
    assert names == sorted(bt.feng_dirs.feng_directions(genes[:7] + ["Absent-syn"])[0]) == sorted(FENG_GROUPS)
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: UNION)
    assert bt._t1_union(True) is UNION
    _feng_stats(stats, genes, dict(cells, blood=30))  # Blood now has a direction
    assert "Blood" in bt._feng_group_names()
    with pytest.raises(SystemExit, match="Blood is among Feng's direction groups"):
        bt._t1_union(True)


def test_build_record_hides_counts_on_a_test_board(bt):
    b = _board()
    rec = {"board": b.key, "n_cells": SECRET_MAX, "record": {"n": SECRET_MAX, "n_before": 2000, "n_after": 2321,
                                                              "n_jittered": 3, "t": 0.5}}
    out = bt._hide_counts(rec, b)
    assert out["n_cells"] == "max_cells" and out["record"] == {"t": 0.5}
    assert str(SECRET_MAX) not in json.dumps(out)
    with pytest.raises(SystemExit, match="does not hold the board's max_cells"):
        bt._hide_counts(dict(rec, n_cells=7), b)
    val = _board(key="T2:heart:val_extrap", split="val_extrap")
    assert bt._hide_counts(rec, val) is rec


def test_build_record_strips_t1s_realised_composition_on_a_test_board(bt):
    """docs/PLAN.md decision log 2026-11-06 (g): build_t1's realised composition is quota / n_cells per category, so
    its denominators give the closed test count; on a test board the record leaves it out (the targets stay)."""
    quotas = {"blood": 1037, "cm": 1146, "endocardium": 2138}  # sum to SECRET_MAX
    rec = {"entry": "M1", "board": "T1:test", "n_cells": SECRET_MAX, "shares": {"blood": 0.24, "cm": 0.26,
                                                                               "endocardium": 0.5},
           "realised": {c: q / SECRET_MAX for c, q in quotas.items()}}
    denominators = [Fraction(v).limit_denominator(10 ** 6).denominator for v in rec["realised"].values()]
    assert math.lcm(*denominators) == SECRET_MAX  # what the stored fractions gave away
    test = Board(key="T1:test", task="T1", genes=GENES, needs_coords=False, split="test", min_cells=10,
                 max_cells=SECRET_MAX)
    out = bt._hide_counts(rec, test)
    assert "realised" not in out and out["n_cells"] == "max_cells" and out["shares"] == rec["shares"]
    assert "realised" in rec  # the builder's own record is not edited in place
    val = Board(key="T1:val", task="T1", genes=GENES, needs_coords=False, split="val", min_cells=10, max_cells=5000)
    assert bt._hide_counts(dict(rec, board="T1:val"), val)["realised"] == rec["realised"]  # validation boards keep it
