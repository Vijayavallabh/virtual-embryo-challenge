"""vec.world_blur (row R69's file WBH and row H37's arms) and its file builder experiments/builds/world_blur.py, on
synthetic arrays and a synthetic h5ad only: nothing is read from data/, no ledger but a temporary one, no score.

Count recovery is exact (smallest counts 1, 2 and 3) and refuses otherwise; a zero spread gives the input's X bytes;
blocks and cell pairs are antithetic; each gene's mean X is held within 1e-3 (fix 2) and the re-measure before the hold
keeps the mean library within 2%; the same bytes twice; the prior flags unexpressed genes and gives them the median
spread (fix 3) and shrinks log s halfway; G0's prior checks, G00 (fix 4a) and the output checks (fix 4b); and every
refusal of the file's builder, then a build that changes X alone."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from vec import world_blur as wb

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "builds"))
import world_blur as cli  # noqa: E402

GENES = [f"g{i}" for i in range(60)]


def counts_and_x(n=900, g=60, seed=3, depth=60.0):
    rng = np.random.default_rng(seed)
    gm = rng.lognormal(-1.0, 1.2, size=g)
    gm *= depth / gm.sum()
    prog = rng.integers(0, 4, size=n)
    lam = gm[None, :] * rng.lognormal(0, 0.7, size=(4, g))[prog]
    c = rng.poisson(lam).astype(float)
    c[:, 0] += 1  # a count-1 gene somewhere in every cell, as the AUDIT found
    x = np.log1p(1e4 * c / c.sum(1, keepdims=True)).astype(np.float32)
    coords = rng.normal(size=(n, 3)) * 20
    return c, x, coords


# ------------------------------------------------------------------------------------------------ counts
def test_count_recovery_is_exact_and_tries_two_and_three():
    c, x, _ = counts_and_x(n=50)
    got, N, used = wb.recover_counts(x)
    assert np.array_equal(got, c) and np.allclose(N, c.sum(1)) and set(used) == {1}
    c2 = c.copy()
    c2[:5] = 0
    c2[:5, :3] = [2, 3, 7]  # smallest count 2: a whole ratio only at 2
    c2[5:8] = 0
    c2[5:8, :3] = [3, 4, 5]  # smallest count 3
    x2 = np.log1p(1e4 * c2 / c2.sum(1, keepdims=True)).astype(np.float32)
    got, N, used = wb.recover_counts(x2)
    assert np.array_equal(got, c2) and list(used[:8]) == [2] * 5 + [3] * 3


def test_count_recovery_refuses_a_cell_that_is_not_whole():
    _, x, _ = counts_and_x(n=20)
    bad = x.copy()
    bad[3, 1] = bad[3, 1] + 0.37
    with pytest.raises(wb.CountRecoveryError):
        wb.recover_counts(bad)
    with pytest.raises(wb.CountRecoveryError):
        wb.recover_counts(-x)


def test_an_empty_cell_recovers_as_empty():
    _, x, _ = counts_and_x(n=20)
    x[4] = 0
    got, N, used = wb.recover_counts(x)
    assert got[4].sum() == 0 and used[4] == 0


# ------------------------------------------------------------------------------------------------ worlds
def test_blocks_are_paired_by_size_and_antithetic():
    _, _, coords = counts_and_x(n=400)
    w = wb.block_worlds(coords, 7, 60)
    assert len(w.pairs) == 16 and w.z.shape == (16, 60) and w.antithetic()
    assert set(np.unique(w.unit)) == set(range(16)) and set(np.unique(w.sign)) == {-1.0, 1.0}
    s = np.full(60, 0.4)
    F = w.factors(s)
    for j in range(16):
        a, b = F[(w.unit == j) & (w.sign > 0)][0], F[(w.unit == j) & (w.sign < 0)][0]
        assert np.allclose(a * b, np.exp(-s ** 2))  # +z and -z: the product is exp(-s^2)
    w.pairs[0] = (w.pairs[0][0], w.pairs[0][0])
    assert not w.antithetic()


def test_cell_pairs_leave_an_odd_cell_out_unchanged():
    w = wb.cell_worlds(9, 5, 10)
    assert len(w.z) == 4 and (w.unit == -1).sum() == 1 and w.antithetic()
    F = w.factors(np.full(10, 0.5))
    assert np.all(F[w.unit == -1] == 1)


def test_the_permuted_spread_is_a_permutation():
    s = np.linspace(0.1, 0.9, 60)
    p = wb.permuted(s, 11)
    assert np.array_equal(np.sort(p), np.sort(s)) and not np.array_equal(p, s)
    assert np.array_equal(p, wb.permuted(s, 11))


def test_too_few_cells_for_the_blocks_are_refused():
    with pytest.raises(ValueError):
        wb.block_worlds(np.zeros((10, 3)), 1, 5)


# ------------------------------------------------------------------------------------------------ the blur
def test_a_zero_spread_gives_the_input_bytes_dense_and_sparse():
    _, x, coords = counts_and_x()
    out = wb.blur(x, coords, np.zeros(60), 1000001)
    assert out.X.dtype == x.dtype and out.X.tobytes() == x.tobytes() and out.record["changed_cells"] == 0
    xs = sparse.csr_matrix(x)
    outs = wb.blur(xs, coords, np.zeros(60), 1000001)
    assert sparse.issparse(outs.X) and (outs.X != xs).nnz == 0


@pytest.fixture(scope="module")
def blurred():
    _, x, coords = counts_and_x()
    s = np.full(60, 0.3)
    return x, coords, s, wb.blur(x, coords, s, 1000001)


def test_the_blur_holds_each_gene_and_the_design_depth(blurred):
    x, coords, s, out = blurred
    r = out.record
    assert r["recovery"]["exact"] and r["antithetic"] and r["worlds"] == 16
    assert np.abs(wb.pseudobulk(out.X) - wb.pseudobulk(x)).max() <= wb.PB_TOL
    assert abs(r["library_ratio_design"] - 1) <= wb.LIBRARY_TOL
    assert r["changed_cells"] >= 0.5 and sum(r["genes"].values()) == 60
    got, _, _ = wb.recover_counts(out.X)  # the output is whole counts over each cell's fixed N
    assert np.all(got >= 0)


def test_unchanged_entries_keep_their_bytes(blurred):
    x, coords, s, out = blurred
    c, N, _ = wb.recover_counts(x)
    c2, _, _ = wb.recover_counts(out.X)
    same = c2 == c
    assert np.array_equal(out.X[same], x[same])


def test_the_same_bytes_twice_and_another_seed_differs(blurred):
    x, coords, s, out = blurred
    again = wb.blur(x, coords, s, 1000001)
    assert again.X.tobytes() == out.X.tobytes()
    assert wb.blur(x, coords, s, 1000002).X.tobytes() != out.X.tobytes()


def test_the_three_arms_differ(blurred):
    x, coords, s, out = blurred
    s2 = np.linspace(0.1, 0.6, 60)
    w = wb.blur(x, coords, s2, 1000001)
    u = wb.blur(x, coords, s2, 1000001, mode=wb.CELLS)
    p = wb.blur(x, coords, s2, 1000001, permute=True)
    assert len({w.X.tobytes(), u.X.tobytes(), p.X.tobytes()}) == 3
    assert u.record["mode"] == wb.CELLS and p.record["permuted"]
    with pytest.raises(ValueError):
        wb.blur(x, coords, s2, 1, mode="other")
    with pytest.raises(ValueError):
        wb.blur(x, coords, -s2, 1)


def test_the_hold_takes_the_plateau_nearest_the_input():
    rng = np.random.default_rng(0)
    n = 4000
    c = rng.poisson(1.0, n).astype(float)
    N = np.full(n, 60.0)
    fp = np.where(np.arange(n) % 2 == 0, 0.7, 1.3)
    lam = np.full(n, 1.0)
    out, how, design = wb.hold_gene(c, fp, lam, N, np.ones(n, dtype=bool), np.random.default_rng(1))
    assert how in ("topup", "thin")
    assert abs(np.log1p(1e4 * out / N).mean() - np.log1p(1e4 * c / N).mean()) <= wb.PB_TOL
    frozen = np.zeros(n, dtype=bool)
    out, how, _ = wb.hold_gene(c, fp, lam, N, frozen, np.random.default_rng(1))
    assert how == "topup" and np.array_equal(out, c)  # no live cell: nothing moves


# ------------------------------------------------------------------------------------------------ the prior
def _table(n, means):
    return {g: (n, np.asarray(m, dtype=float)) for g, m in means.items()}


def test_the_prior_flags_unexpressed_genes_and_shrinks_log_s():
    ng = 6
    base = np.array([5.0, 8.0, 0.2, 0.1, 12.0, 3.0])
    rng = np.random.default_rng(0)
    tables = {}
    for i, f in enumerate(wb.PRIOR_FILES):
        tables[f] = _table(300, {"NCC": base * np.exp(rng.normal(0, 0.3, ng)) * (1 + 0.1 * i),
                                 "Neural": base * np.exp(rng.normal(0, 0.3, ng)),
                                 "Endoderm": base * np.exp(rng.normal(0, 0.3, ng))})
        tables[f]["Ectoderm"] = (100, base)  # under 200 cells: never used
    for f in tables:
        for g in ("NCC", "Neural", "Endoderm"):
            tables[f][g][1][2:4] = 0.3  # genes 2 and 3 under 1 CP10k everywhere: flagged
    p = wb.prior_from_tables(tables, [f"g{i}" for i in range(ng)])
    assert p.used == list(wb.PAIR_NAMES) and all(len(v) == 3 for v in p.groups.values())
    assert list(p.flagged) == [False, False, True, True, False, False]
    keep = ~p.flagged
    med = np.median(np.log(np.maximum(p.raw[keep], 1e-6)))
    assert np.allclose(np.log(p.s[keep]), 0.5 * np.log(p.raw[keep]) + 0.5 * med)
    assert np.allclose(p.s[p.flagged], np.exp(med))
    sm = p.summary()
    assert sm["flagged"] == 2 and sm["n_genes"] == ng and set(sm) >= {"median_s", "ko_wt_median_spread"}


def test_a_pair_with_one_group_is_not_used_and_no_pair_gives_nan():
    one = {f: _table(300, {"NCC": np.ones(3)}) for f in wb.PRIOR_FILES}
    p = wb.prior_from_tables(one, ["a", "b", "c"])
    assert p.used == [] and np.isnan(p.s).all()
    assert wb.prior_checks(p) == {"pairs": False, "median s": False, "ko against replicate": False}


def test_the_prior_checks():
    g = ("a", "b")
    prof = {n: np.array([0.3, -0.3]) for n in wb.PAIR_NAMES}
    p = wb.Prior(g, np.array([0.3, 0.3]), np.array([0.3, 0.3]), np.zeros(2, dtype=bool), prof, {})
    assert all(wb.prior_checks(p).values())
    p.s = np.array([1.5, 1.5])
    assert not wb.prior_checks(p)["median s"]
    big = {**prof, **{n: np.array([2.0, -2.0]) for n in wb.KO_WT}}
    p2 = wb.Prior(g, np.array([0.3, 0.3]), np.array([0.3, 0.3]), np.zeros(2, dtype=bool), big, {})
    assert not wb.prior_checks(p2)["ko against replicate"]
    three = {n: prof[n] for n in wb.PAIR_NAMES[:2]}
    p3 = wb.Prior(g, np.array([0.3, 0.3]), np.array([0.3, 0.3]), np.zeros(2, dtype=bool), three, {})
    assert not wb.prior_checks(p3)["pairs"]


def test_group_means_read_non_mesp1_tracked_groups_only():
    C = sparse.csr_matrix(np.array([[1, 1], [2, 0], [0, 0], [3, 1]], dtype=float))
    out = wb.group_means(C, ["NCC", "NCC", "NCC", "Untracked"], ["no", "no", "no", "no"])
    assert list(out) == ["NCC"] and out["NCC"][0] == 2  # the empty cell and the untracked group left out
    assert np.allclose(out["NCC"][1], [7500, 2500])
    assert wb.group_means(C, ["CM"] * 4, ["yes"] * 4) == {}


# ------------------------------------------------------------------------------------------------ G00 and checks
def test_g00_passes_a_wide_world_and_stops_a_narrow_one():
    _, x, _ = counts_and_x(n=500)
    wide = wb.g00(x, np.full(60, 3.0), 1000001)
    narrow = wb.g00(x, np.full(60, 0.01), 1000001)
    assert wide["pass"] and wide["cut"] >= wb.G00_MIN_CUT and not narrow["pass"]
    assert narrow["c"] < wide["c"] and wb.g00(x, np.full(60, 3.0), 1000001) == wide
    assert not wb.g00(x, np.full(60, 1.0), 1000001)["pass"]  # sparse cells under a wide median kernel


def test_the_output_checks(blurred):
    x, coords, s, out = blurred
    refs = {"r": wb.pseudobulk(x) * 0.9}
    assert all(wb.output_checks(x, out, refs).values())
    off = wb.Blurred(out.X, {**out.record, "library_ratio_design": 1.05, "antithetic": False})
    ch = wb.output_checks(x, off, refs)
    assert not ch["library"] and not ch["antithetic"] and ch["count recovery"]
    assert not wb.output_checks(x, wb.Blurred(x, {**out.record, "changed_cells": 0.0}), refs)["changed cells"]
    gone = wb.output_checks(x, None, refs, "CountRecoveryError")
    assert not any(gone.values())
    held = wb.output_checks(x, None, refs, "HoldError")
    assert held["count recovery"] and not held["pseudobulk rank"]


# ------------------------------------------------------------------------------------------------ the file's builder
R69 = ("R69 (the captain's order of 2026-11-25, ruling, no run: \"File WBH may use ... as inputs to its world prior "
       "only\", a synthetic stand-in; built only after its gate row H37 reads PASS, only with `--captain-ruling R69`;"
       " not a try)")
FIRST = ("H37 (the gate of row R69's file WBH, the heart's world-blur hedge: six counted tries, arms W, U and P "
         "on two transitions; G0 first, target-free, not a try)")


def ledger(tmp_path, outcome="OUTCOME H37: PASS (world pattern).", r69=R69, first=FIRST) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text(f"| {r69} | heart | x |\n| {first} | heart | x | x | x | x | x | x | x | x | {outcome} |\n")
    return p


def half_file(tmp_path, name=cli.IN_NAME, n=400) -> Path:
    _, x, coords = counts_and_x(n=n)
    a = ad.AnnData(X=x, var={"g": GENES}, obs={"o": np.arange(n)})
    a.var_names = GENES
    a.obs_names = [f"cell_{i}" for i in range(n)]
    a.obsm["spatial_3D"] = coords.astype(np.float32)
    d = tmp_path / "in"
    d.mkdir(exist_ok=True)
    p = d / name
    a.write_h5ad(p, compression="gzip")
    return p


def sums(tmp_path, path: Path, digest: str | None = None) -> Path:
    s = tmp_path / "SUMS"
    d = digest or hashlib.sha256(path.read_bytes()).hexdigest()
    s.write_text(f"{'0' * 64}  tka/T2heart_organ_ts125_tk.h5ad\n{d}  tka/{cli.IN_NAME}\n{d}  tkb/{cli.IN_NAME}\n")
    return s


def argv(src, out, **over):
    a = {"--prior": "same_stage", "--blocks": "32", "--seed": "1000001", "--captain-ruling": "R69", **over}
    return [str(src), str(out), *[x for k, v in a.items() for x in (k, v)]]


def _refused(tmp_path, args, tries, sm) -> str:
    with pytest.raises(SystemExit) as e:
        cli.refusals(cli.build_parser().parse_args(args), tries, sm)
    return str(e.value)


def test_the_builder_constants_and_record():
    assert cli.RULING == "R69" and cli.SEED == 1000001 and cli.BLOCKS == 32 and cli.PRIOR == "same_stage"
    assert cli.recorded_sha256() == "9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54"


def test_every_refusal_of_the_builder(tmp_path):
    src = half_file(tmp_path)
    sm, tries = sums(tmp_path, src), ledger(tmp_path)
    out = tmp_path / "o" / "wb.h5ad"
    cli.refusals(cli.build_parser().parse_args(argv(src, out)), tries, sm)  # nothing refused
    assert "--captain-ruling R69" in _refused(tmp_path, argv(src, out, **{"--captain-ruling": "R55"}), tries, sm)
    for k, v in (("--prior", "other"), ("--blocks", "16"), ("--seed", "1000002")):
        assert "constants are fixed" in _refused(tmp_path, argv(src, out, **{k: v}), tries, sm)
    no_r69 = tmp_path / "n1"
    no_r69.mkdir()
    assert "row R69" in _refused(tmp_path, argv(src, out), ledger(no_r69, r69="R69 (other) x"), sm)
    for word in ("OUTCOME H37: STOP: IDLE.", "OUTCOME H37: VOID: x; nothing is built.", "pending",
                 "OUTCOME H37: PASS (world pattern). OUTCOME H37: STOP: IDLE."):
        d = tmp_path / f"w{abs(hash(word))}"
        d.mkdir()
        assert "PASS line" in _refused(tmp_path, argv(src, out), ledger(d, outcome=word), sm)
    d = tmp_path / "first"
    d.mkdir()
    assert "PASS line" in _refused(tmp_path, argv(src, out), ledger(d, first="H37 (another first cell)"), sm)
    assert "IN must be" in _refused(tmp_path, argv(src, out), tries, sums(tmp_path, src, "f" * 64))
    (tmp_path / "x").mkdir()
    other = half_file(tmp_path / "x", name="other.h5ad")
    assert "IN must be" in _refused(tmp_path, argv(other, out), tries, sums(tmp_path, src))
    assert "OUT must not be IN" in _refused(tmp_path, argv(src, src), tries, sums(tmp_path, src))
    out.parent.mkdir()
    out.write_text("x")
    assert "exists" in _refused(tmp_path, argv(src, out), tries, sums(tmp_path, src))


def test_the_builder_changes_x_alone_and_gives_the_same_bytes_twice(tmp_path, monkeypatch, own_repo):
    src = half_file(tmp_path)
    sm, tries = sums(tmp_path, src), ledger(tmp_path)
    prior = wb.Prior(tuple(GENES), np.full(60, 0.3), np.full(60, 0.3), np.zeros(60, dtype=bool), {}, {})
    monkeypatch.setattr(cli.wb, "load_prior", lambda data, genes: prior)
    outs = [tmp_path / "a" / "wb.h5ad", tmp_path / "b" / "wb.h5ad"]
    digests = [cli.main(argv(src, o), tries, sm, GENES) for o in outs]
    assert digests[0] == digests[1] and outs[0].read_bytes() == outs[1].read_bytes()
    a, b = ad.read_h5ad(src), ad.read_h5ad(outs[0])
    assert list(a.obs_names) == list(b.obs_names) and list(a.var_names) == list(b.var_names)
    assert a.obsm["spatial_3D"].tobytes() == b.obsm["spatial_3D"].tobytes() and b.X.dtype == a.X.dtype
    assert not np.array_equal(a.X, b.X)
    assert (outs[0].parent / "wb.h5ad.stamp.json").exists()
    with pytest.raises(SystemExit):
        cli.main(argv(src, tmp_path / "c" / "wb.h5ad"), tries, sm, list(reversed(GENES)))


def test_the_ruling_check_on_the_committed_ledger():
    """Row R69 as committed is the captain's order the builder needs (read by code; nothing printed)."""
    rows = cli._rows()
    assert rows.ruling_ok((REPO / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8"))


def test_refusal_args_namespace_shape():
    a = SimpleNamespace(captain_ruling="R69", prior="same_stage", blocks=32, seed=1000001, src="x", out="y")
    with pytest.raises(SystemExit):
        cli.refusals(a, Path("/nonexistent/TRIES.md"), Path("/nonexistent/SUMS"))
