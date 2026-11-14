"""T1 arm `b0_xt` (campaign/t1-d): b0 on T1 with its momentum step, B0's only count edit on T1, realised by the
exact-total kernel in its repair form (heart arm b0_trend_xt's kernel, `vec.core.apply_factors(kernel="exact_total")`,
unchanged), with its check line and record. Synthetic data only; nothing scored, no prediction compared with a
target."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import core, harness, t1_scope
from vec.boards import Board
from vec.harness import B0_DAMPING, CANDIDATES, BuildContext, PseudoTask, TaskView, load_tasks
from vec.prep import log1p_cp10k

GENES = [f"g{i}" for i in range(50)]
LABELS = ["V-CM", "Endocardium", "Blood", "Macrophage"]
HELD = ("g49",)
TASK = PseudoTask(key="T1/E9.5->E10.5", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                  reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                  validation_board="T1:val")


def _labelled(n, seed):
    rng = np.random.default_rng(seed)
    mu = rng.gamma(2.0, 8.0, len(GENES))
    C = rng.poisson(mu, (n, len(GENES))).astype(np.int32)
    C[np.arange(n), rng.integers(0, len(GENES), n)] += 1
    C[: n // 3, 5] = 0  # a gene with zeros in many cells
    a = ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), layers={"counts": sparse.csr_matrix(C)},
                   obs=pd.DataFrame({"celltype": np.resize(LABELS, n)}, index=[f"c{i}" for i in range(n)]))
    a.var_names = GENES
    return a


@pytest.fixture(scope="module")
def stages():
    return {"ref": _labelled(1500, 1), "start": _labelled(1200, 2)}


def _ctx(inputs, seed=harness.BUILD_SEED_OFFSET + 1):
    board = Board(key="pseudo", task="T1", genes=list(GENES), needs_coords=False, min_cells=50, max_cells=700)
    return BuildContext(task=TaskView.of(TASK), board=harness._candidate_board(board), inputs=inputs, seed=seed,
                        derived={"t1_held_flat": HELD})


def _spy(monkeypatch):
    calls = []
    real = t1_scope.apply_factors

    def spy(C, F, tidx, held=None, keep_depth=True, seed=0, **kw):
        out = real(C, F, tidx, held, keep_depth=keep_depth, seed=seed, **kw)
        calls.append({"C": C.copy(), "F": F.copy(), "tidx": tidx.copy(), "seed": list(np.atleast_1d(seed)),
                      "kernel": kw.get("kernel", "independent"), "keep_depth": keep_depth, "out": out[0].copy(),
                      "scale": out[1].copy(), "dead": out[2].copy(),
                      "held": np.zeros(C.shape[1], bool) if held is None else np.asarray(held, bool).copy()})
        return out

    monkeypatch.setattr(t1_scope, "apply_factors", spy)
    return calls


def _flat(stages):
    """start = the reference's own cells: every type's momentum is exactly 0, every factor exactly 1."""
    return {"ref": stages["ref"], "start": stages["ref"].copy()}


# ------------------------------------------------------------------------------------------------- 1. b0 is unchanged
def test_kernel_unset_b0_is_unchanged(stages, monkeypatch):
    """b0 makes one apply_factors call on T1 (the momentum step), with the default kernel and no repair record; the
    branch's pins of b0 and of the shared core run in tests/test_b0_trend_mix_xt.py and tests/test_core.py."""
    calls = _spy(monkeypatch)
    CANDIDATES["b0"].build(_ctx(stages))
    assert len(calls) == 1 and calls[0]["kernel"] == "independent"


def test_repair_record_records_and_changes_nothing():
    rng = np.random.default_rng(31)
    C = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 2.0, 40), (60, 40)))
    F = np.exp(rng.normal(0, 0.6, (3, 40)))
    ti = rng.integers(0, 3, 60)
    rec: dict = {}
    a = core.apply_factors(C, F, ti, seed=[5, 3], kernel="exact_total")
    b = core.apply_factors(C, F, ti, seed=[5, 3], kernel="exact_total", repair_record=rec)
    for x, y in zip(a, b, strict=True):
        dense_x, dense_y = (z.toarray() if sparse.issparse(z) else z for z in (x, y))
        np.testing.assert_array_equal(dense_x, dense_y)
    y_ind = core.apply_factors(C, F, ti, seed=[5, 3])[0].toarray()
    D = C.toarray().astype(np.float64)
    f = F[ti] * a[1][:, None]
    M = (f != 1.0) & (D > 0)
    np.testing.assert_array_equal(rec["off"], (y_ind * M).sum(1) - np.rint((D * f * M).sum(1)).astype(np.int64))
    np.testing.assert_array_equal(rec["moved_cells"], M.any(1))
    empty: dict = {}
    core.apply_factors(C, F, ti, seed=[5, 3], repair_record=empty)
    assert empty == {}  # the independent kernel records nothing


# ------------------------------------------------------------------------------------ 2. the kernel on a moving step
def test_b0_xt_is_b0s_step_with_each_moved_total_at_its_rounded_expected_total(stages, monkeypatch):
    calls = _spy(monkeypatch)
    ctx = _ctx(stages)
    CANDIDATES["b0"].build(ctx)
    CANDIDATES["b0_xt"].build(ctx)
    assert len(calls) == 2
    b0, xt = calls
    np.testing.assert_array_equal(xt["C"].toarray(), b0["C"].toarray())
    np.testing.assert_array_equal(xt["F"], b0["F"])
    np.testing.assert_array_equal(xt["tidx"], b0["tidx"])
    np.testing.assert_array_equal(xt["held"], b0["held"])
    np.testing.assert_array_equal(xt["scale"], b0["scale"])
    assert xt["seed"] == b0["seed"] and xt["keep_depth"] and b0["kernel"] == "independent"
    assert xt["kernel"] == "exact_total"
    assert not xt["dead"].any() and not b0["dead"].any()
    D, y, x = b0["C"].toarray().astype(np.float64), b0["out"].toarray(), xt["out"].toarray()
    f = np.where(b0["held"][None, :], 1.0, b0["F"][b0["tidx"]] * b0["scale"][:, None])
    M = (f != 1.0) & (D > 0)
    assert (f != 1.0).any(), "the synthetic momentum is not flat"
    target = np.rint((D * f * M).sum(1)).astype(np.int64)
    d = (y * M).sum(1).astype(np.int64) - target
    np.testing.assert_array_equal((x * M).sum(1), target)  # every moved cell's total at its rounded expected total
    np.testing.assert_array_equal(x[~M], y[~M])  # entries the step does not move: as b0 leaves them
    held = np.isin(GENES, HELD)
    assert held.any() and np.array_equal(x[:, held], D[:, held])  # held genes untouched
    assert (x[D == 0] == 0).all()  # zeros stay zero
    assert (d != 0).any() and (d == 0).any() and np.array_equal(x[d == 0], y[d == 0])
    np.testing.assert_array_equal(np.abs(x - y).sum(1), np.abs(d))
    rec = harness.T1_ARM_RECORDS["b0_xt"]
    moved = M.any(1)
    assert rec["repaired_share"] == pytest.approx(float((d[moved] != 0).mean()), abs=0, rel=0)
    assert rec["mean_abs_off"] == pytest.approx(float(np.abs(d[moved]).mean()), abs=0, rel=0)


# ----------------------------------------------------------------------------------------------- 3. a flat momentum
def test_flat_momentum_is_b0_bit_for_bit(stages, monkeypatch):
    calls = _spy(monkeypatch)
    inputs = _flat(stages)
    a, b = CANDIDATES["b0_xt"].build(_ctx(inputs)), CANDIDATES["b0"].build(_ctx(inputs))
    xt = calls[0]
    assert (xt["F"] == 1.0).all() and (xt["scale"] == 1.0).all() and xt["kernel"] == "exact_total"
    np.testing.assert_array_equal(a.X, b.X)
    rec = harness.T1_ARM_RECORDS["b0_xt"]
    assert rec["repaired_share"] == 0.0 and rec["mean_abs_off"] == 0.0


# ------------------------------------------------------------------------------------------- 4. the realised library
LIBRARY_RTOL = 0.005  # the test's tolerance on the mean realised library, b0_xt against b0 (not a constant of the arm)


def test_mean_realised_library_is_b0s(stages, monkeypatch):
    """keep_depth makes each cell's expected total its original total; b0_xt realises it exactly per cell (the
    rounded expected total of a sum that is an integer), b0 on average; the means agree within LIBRARY_RTOL."""
    calls = _spy(monkeypatch)
    ctx = _ctx(stages)
    CANDIDATES["b0"].build(ctx)
    CANDIDATES["b0_xt"].build(ctx)
    b0, xt = calls
    lib_in = b0["C"].toarray().sum(1)
    lib_b0, lib_xt = b0["out"].toarray().sum(1), xt["out"].toarray().sum(1)
    np.testing.assert_array_equal(lib_xt, lib_in)
    assert lib_b0.mean() == pytest.approx(lib_xt.mean(), rel=LIBRARY_RTOL)
    assert not np.array_equal(lib_b0, lib_in)


# --------------------------------------------------------------------------------------------------------- 5. seeds
def test_same_seed_same_output_another_seed_another_draw(stages):
    a = CANDIDATES["b0_xt"].build(_ctx(stages))
    b = CANDIDATES["b0_xt"].build(_ctx(stages))
    c = CANDIDATES["b0_xt"].build(_ctx(stages, seed=harness.BUILD_SEED_OFFSET + 2))
    np.testing.assert_array_equal(a.X, b.X)
    assert a.X.shape == c.X.shape and not np.array_equal(a.X, c.X)


# ---------------------------------------------------------------------------------------------------- 6. check line
CHECK = re.compile(r"b0_xt check: task=(\S+) seed=(\d+) repaired_share=(\d+\.\d{4}) mean_abs_off=(\d+\.\d{4})")


def test_check_line_holds_only_the_allowed_fields(stages, capsys):
    ctx = _ctx(stages)
    CANDIDATES["b0_xt"].build(ctx)
    lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.startswith("b0_xt check:")]
    assert len(lines) == 1
    m = CHECK.fullmatch(lines[0])
    assert m and m.group(1) == TASK.key and int(m.group(2)) == ctx.seed
    rec = harness.T1_ARM_RECORDS["b0_xt"]
    assert set(rec) == {"arm", "arm_rule", "kernel", "task", "seed", "repaired_share", "mean_abs_off"}
    assert float(m.group(3)) == pytest.approx(rec["repaired_share"], abs=5e-5)
    assert float(m.group(4)) == pytest.approx(rec["mean_abs_off"], abs=5e-5)
    assert 0.0 < rec["repaired_share"] <= 1.0 and rec["mean_abs_off"] > 0.0
    assert not any(g in lines[0] for g in GENES)


# ---------------------------------------------------------------------------------------------------- registration
def test_registered_for_t1_with_b0s_tasks_and_requirements():
    tasks = [t for ts_ in load_tasks().values() for t in ts_]
    c, base = CANDIDATES["b0_xt"], CANDIDATES["b0"]
    assert c.tasks == frozenset({"T1"}) and c.requires is base.requires and c.derived == base.derived
    assert c.needs == base.needs
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if base.applies(t) and t.task == "T1"]
    assert any(c.applies(t) for t in tasks)
    assert c.params["scope"] == "chamber" and c.params["damping"] == B0_DAMPING
    assert c.params["rule"] == "vec.harness.B0_XT_RULE" and isinstance(harness.B0_XT_RULE, str)


def test_the_arm_refuses_other_combinations():
    with pytest.raises(ValueError, match="unknown kernel"):
        harness._b0_t1(B0_DAMPING, kernel="multinomial")
    with pytest.raises(ValueError, match="B0 alone"):
        harness._b0_t1(0.0, kernel="exact_total")
    with pytest.raises(ValueError, match="B0 alone"):
        harness._b0_t1(B0_DAMPING, herd=True, kernel="exact_total")
    with pytest.raises(ValueError, match="B0 alone"):
        harness._b0_t1(B0_DAMPING, prior=lambda ctx: None, kernel="exact_total")
