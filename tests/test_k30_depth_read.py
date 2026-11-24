"""Try K30 (experiments/harness/TRIES.md; experiments/t3_decomp/k30_depth_read.py): the read to which ruling row R53
binds one T3 upload, the frozen per-label depth step on K18's barred_form. Synthetic vocabularies, pools and depth
tables, hand-made per-seed score tables, hand-made task files and hand-made K18 files whose pins are computed here, a
hand-made ledger and a fixed clock. No competition file, no score of a fold, no read of the committed ledger and no
value of the committed depth table (only its bytes are hashed against the script's pin). The launcher's tests are in
tests/test_run_k30.py."""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3 = REPO / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import k7_wt_revert as k7  # noqa: E402
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402
import k30_depth_read as k30  # noqa: E402

SCEN = ["OFT/RV-CM", "aSHF", "pSHF", "aPHM", "pPHM"]
LABELS = ["vCM1", "aCM1", "OFT/RV-CM", "aSHF", "pSHF", "aPHM", "pPHM", "Endo", "Blood"]
SPEC = {"factor": {"source": "mab21l2", "exponent": 1.0},
        "composition": {"kind": "central", "depleted_states": SCEN, "spared_states": ["vCM1"],
                        "log2_magnitude": 0.717813870986706, "depletion_factor": 2.0 ** -0.717813870986706,
                        "magnitude_source": "Mab21l2 (D7m)"},
        "depth": False, "stage": False, "coverage": False}
# The synthetic depth table: rep1 and rep2 ratios and WT cells per label. vCM1 is at 1.0 in rep2 (not listed for g1),
# aCM1 at 1.0 in rep1 (not listed for g2), aSHF has no WT cells (weight 0), and Gone, at 1.0 in both columns, is not a
# label of the test vocabularies (a table label the vocabulary lacks is harmless when it is not listed)
TABLE_ROWS = [  # label, wt_cells, rep1_depth_ratio, rep2_depth_ratio
    ("vCM1", 100, 0.5, 1.0), ("aCM1", 50, 1.0, 0.4), ("OFT/RV-CM", 25, 2.0, 0.8), ("aSHF", 0, 0.25, 1.25),
    ("pSHF", 40, 0.6, 0.7), ("aPHM", 10, 0.9, 0.5), ("pPHM", 10, 1.5, 0.6), ("Endo", 30, 0.7, 0.75),
    ("Blood", 20, 0.8, 0.9), ("Gone", 15, 1.0, 1.0)]
SCRIPT_P, LAUNCHER_P = "`experiments/t3_decomp/k30_depth_read.py`", "`experiments/t3_decomp/run_k30.sh"
ROW = ("| K30 (one row, one script `experiments/t3_decomp/k30_depth_read.py` and its launcher "
       "`experiments/t3_decomp/run_k30.sh <COMMIT> <RUN>`; S3 weight `--s3-weight plain`; three tries) | T3 | the "
       "row's other cells |\n")
R53 = ("| R53 (the captain's order of 2026-11-24, ruling, no run: \"Yes: ...\"; row K30's world weight is `plain`; not "
       "a try) | T3 | the ruling's other cells |\n")
LEDGER = "| K29 (another row) | T3 | x |\n" + R53 + ROW
G, N_POOL = 12, 120
GENES = [f"g{i}" for i in range(G)]
PLAIN, DM = k30.PLAIN_CELLS, k30.DM_CELLS
BEFORE = datetime(2026, 11, 25, 6, 0, tzinfo=UTC)
RESULTS_AS_SHIPPED = k30.RESULTS


def _table_text(rows=TABLE_ROWS, extra_col=True) -> str:
    head = "label,group,wt_cells,wt_share,rep1_depth_ratio,rep2_depth_ratio" if extra_col else \
        "label,wt_cells,rep1_depth_ratio,rep2_depth_ratio"
    body = [f"{t},grp,{w},0.1,{r1!r},{r2!r}" if extra_col else f"{t},{w},{r1!r},{r2!r}" for t, w, r1, r2 in rows]
    return "\n".join([head, *body]) + "\n"


@pytest.fixture(autouse=True)
def synthetic_table(monkeypatch, tmp_path):
    """Every test reads a synthetic table at its own pin, reads no clock (before the cutoff), meets no reading on disk
    and leaves no task set."""
    p = tmp_path / "depth_table.csv"
    p.write_text(_table_text())
    monkeypatch.setattr(k30, "TABLE", p)
    monkeypatch.setattr(k30, "TABLE_SHA256", hashlib.sha256(p.read_bytes()).hexdigest())
    monkeypatch.setattr(k30, "RESULTS", tmp_path / "k30_results_none")
    monkeypatch.setattr(k30, "_now", lambda: BEFORE)
    yield p
    k30._K30.clear()
    k18._FOLD.clear()


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _vocab(labels=None, seed=0):
    t3d = harness._t3_head().t3d
    labels = list(reversed(LABELS)) + ["Extra"] if labels is None else labels  # another order than the table's
    count = np.random.default_rng(seed).integers(20, 200, len(labels))
    L = len(labels)
    return t3d.Vocab(labels=list(labels), group=np.array(["g"] * L), mesp1=np.array(["partial"] * L),
                     nm=np.zeros(L, dtype=bool), broad=np.ones(L, dtype=bool), model_idx=np.arange(L),
                     p=count / count.sum(), count=count)


def _want(task: str) -> dict[str, float]:
    """The other replicate's ratio per label by hand: g1 takes rep2's column, g2 takes rep1's."""
    col = 3 if task == "g1" else 2
    return {r[0]: r[col] for r in TABLE_ROWS}


# ------------------------------------------------------------------------------------------------ constants, the rows


def test_arms_constants_and_outcome_lines():
    assert k30.ARMS == ("barred_form", "dx1", "dxh", "dxg") and k30.BASE == k18.BARRED == "barred_form"
    assert k30.COUNTED == k30.DEPTH_ARMS == ("dx1", "dxh", "dxg") and k30.EXPONENT == {"dx1": 1.0, "dxh": 0.5}
    assert k30.OTHER == {"g1": "rep2", "g2": "rep1"} and k30.OTHER_TASK == {"g1": "g2", "g2": "g1"}
    assert k30.COLUMN == {"g1": "rep2_depth_ratio", "g2": "rep1_depth_ratio"} and k30.WEIGHTS == "wt_cells"
    assert (k30.NAME_GAIN, k30.MIN_AHEAD, k30.MAX_EXTRA_GATED, k30.NULL_BAND) == (1.0, 8, 1, 0.3)
    assert k30.EV_BAR == 1 and k30.S3_WEIGHTS == {"0.5": Fraction(1, 2), "plain": Fraction(1)}
    assert (k30.GUARD_MEAN, k30.GUARD_SEED) == (0.1, 0.5) and k30.WEIGHT_OPTION == "--s3-weight"
    assert k30.CUTOFF == datetime(2026, 11, 25, 13, 0, tzinfo=UTC) and k30.CUTOFF_TEXT == "2026-11-25 13:00 UTC"
    assert k30.STREAM == "t3_s2_form" and k30.CELLS == ("plain g1", "plain g2", "dm g1", "dm g2")
    assert list(SEEDS) == list(range(10)) and k30.TASKS == k18.TASKS
    assert (k30.RULING, k30.RULING_START) == ("R53", "R53 (the captain's ")
    assert k30.RULING_WEIGHT.format("plain") == "row K30's world weight is `plain`"
    assert k30.FILE_OF == {"dx1": "T3_c2_s2_dd.h5ad", "dxh": "T3_c2_s2_ddh.h5ad", "dxg": "T3_c2_s2_dg.h5ad"}
    assert k30.VOID_LINE.format("x") == "Try K30: VOID: x; no file."
    assert k30.STOP_LINE.format("S1, S4") == "Try K30: STOP (S1, S4): no file; C2 stays."
    assert k30.NAME_LINE.format(arm="dxh") == ("Try K30: NAME dxh: the file may be built twice and uploaded under the "
                                               "captain's bound yes (R53).")
    assert k30.LATE == "INCOMPLETE past its cutoff"
    assert (k30.EXIT_REFUSED, k30.EXIT_INCOMPLETE, k30.EXIT_EXISTS, k30.EXIT_GUARD, k30.EXIT_LATE) == (2, 3, 4, 6, 7)
    assert k30.K18_SHA256 == k21.K18_SHA256 and k30.K18_TASKS == T3 / "results" / "k18_built_files" / "tasks"
    assert RESULTS_AS_SHIPPED == T3 / "results" / "k30_depth_read"
    assert k30.SCRIPT == "experiments/t3_decomp/k30_depth_read.py"
    assert k30.LAUNCHER == "experiments/t3_decomp/run_k30.sh"
    assert k30.TABLE_REL == "experiments/t3_decomp/results/gata4/gata4_composition.csv"


def test_the_committed_table_is_at_the_scripts_pin():
    """Only the pin is checked (the committed table's bytes hashed against the script's constant), never a value. The
    autouse fixture points the module at a synthetic table, so the constant is read from the script's text."""
    pin = re.search(r'^TABLE_SHA256 = "([0-9a-f]{64})"$', (T3 / "k30_depth_read.py").read_text(), flags=re.M).group(1)
    assert _sha(REPO / k30.TABLE_REL) == pin


ROW_CASES = [
    (LEDGER, "plain"),
    (R53 + ROW + f"| R60 (a ruling) | cites {SCRIPT_P} and {LAUNCHER_P}` |\n", "plain"),  # other rows do not block
    (R53.replace("`plain`", "`0.5`") + ROW.replace("`--s3-weight plain`", "`--s3-weight 0.5`"), "0.5"),
    (R53 + ROW.replace("S3 weight `--s3-weight plain`", "`--s3-weight plain` twice `--s3-weight plain`"), "plain"),
    (R53 + ROW * 2, None),  # exactly one row K30
    (R53 + ROW + "| K30 | T3 | a second line opening as row K30 |\n", None),
    (R53 + ROW.replace("| K30 (", "| K30a ("), None),  # K30a is another row
    (R53 + ROW.replace("| K30 (", "| K30 [", 1), None),  # the first cell must start 'K30 ('
    (R53 + ROW.replace(SCRIPT_P, "`experiments/t3_decomp/k29_cd7_read.py`"), None),  # it names this script
    (R53 + ROW.replace(LAUNCHER_P, "`experiments/t3_decomp/run_k29.sh"), None),  # and this launcher
    (R53 + ROW.replace("k30_depth_read.py`", "k30_depth_read.py.bak`"), None),  # each as a whole path
    (R53 + ROW.replace(LAUNCHER_P, "`old/experiments/t3_decomp/run_k30.sh"), None),
    (R53 + f"| K30 (one row; script and launcher in cell 3; `--s3-weight plain`) | {SCRIPT_P} {LAUNCHER_P}` |\n",
     None),
    (R53 + "|  K30 (two spaces) " + ROW[len("| K30 ("):], None),
    (R53 + ROW.replace("S3 weight `--s3-weight plain`; ", ""), None),  # no weight named
    (R53 + ROW.replace("`--s3-weight plain`", "--s3-weight plain"), None),  # not a backticked string
    (R53 + ROW.replace("`--s3-weight plain`", "`--s3-weight 0.7`"), None),  # not an allowed value
    (R53 + ROW.replace("`--s3-weight plain`", "`--s3-weight plain` and `--s3-weight 0.5`"), None),  # both
    (R53 + ROW.replace("`--s3-weight plain`", "`--s3-weight plain` and `--s3-weight 0.7`"), None),
    (R53 + ROW.replace("`--s3-weight plain`", "`--s3-weight plain --s3-weight 0.5`"), None),
    (ROW, None),  # no row R53
    (R53 * 2 + ROW, None),
    (R53.replace("| R53 (", "| R53a (") + ROW, None),
    (R53.replace("R53 (the captain's ", "R53 (a ruling of ") + ROW, None),  # its first cell starts so
    (R53.replace("row K30's world weight is `plain`; ", "") + ROW, None),  # and names row K30's weight
    (R53.replace("`plain`", "`0.5`") + ROW, None),  # the same one
    (R53.replace("is `plain`", "is `plain`, row K30's world weight is `0.5`") + ROW, None),
]


@pytest.mark.parametrize("text_,weight", ROW_CASES)
def test_the_row_rule(tmp_path, capsys, text_, weight):
    p = tmp_path / "TRIES.md"
    p.write_text(text_)
    if weight:
        assert k30.row_present(p) == {"id": "K30", "weight": weight, "ruling": "R53"}
    else:
        with pytest.raises(SystemExit) as e:
            k30.row_present(p)
        assert e.value.code == 2 and "refused: " in capsys.readouterr().err


def test_a_missing_ledger_is_refused(tmp_path):
    with pytest.raises(SystemExit) as e:
        k30.row_present(tmp_path / "missing.md")
    assert e.value.code == 2


# ------------------------------------------------------------------------------------------------ the table


def test_the_table_and_the_other_replicate_pairing():
    t = k30.load_table()
    assert t["labels"] == [r[0] for r in TABLE_ROWS] and t["sha256"] == k30.TABLE_SHA256
    assert t["weights"] == {r[0]: float(r[1]) for r in TABLE_ROWS}
    for task in ("g1", "g2"):  # g1 takes rep2's ratios, g2 takes rep1's: another sample of the same run
        assert t["ratios"][task] == _want(task)
        assert k30.listed(t, task) == {x: r for x, r in _want(task).items() if r != 1.0}
    assert "vCM1" not in k30.listed(t, "g1") and "aCM1" not in k30.listed(t, "g2")
    assert k30.listed(t, "g1")["aCM1"] == 0.4 and k30.listed(t, "g2")["vCM1"] == 0.5


def test_dxg_is_the_wt_cell_weighted_geometric_mean_over_all_labels_one_where_none():
    t = k30.load_table()
    for task in ("g1", "g2"):
        rows = [(w, rr[1] if task == "g1" else rr[0]) for _, w, *rr in TABLE_ROWS]
        num = sum(w * math.log(r) for w, r in rows)  # the 1.0 labels weigh in at log 1 = 0; aSHF's weight is 0
        assert t["global"][task] == pytest.approx(math.exp(num / sum(w for w, _ in rows)), rel=1e-14)
    # over the labelled ratios only it would be another number: the 1.0 labels are in
    only = [(w, r) for _, w, _, r in TABLE_ROWS if r != 1.0]
    assert t["global"]["g1"] != pytest.approx(math.exp(sum(w * math.log(r) for w, r in only) / sum(w for w, _ in only)))
    # by hand: weights 3 and 1, ratios 0.25 and 1.0: exp((3 log 0.25 + 0) / 4) = 0.25^(3/4); unweighted 0.5
    assert k30.weighted_geomean({"a": 0.25, "b": 1.0}, {"a": 3, "b": 1}) == pytest.approx(0.25 ** 0.75, rel=1e-14)
    assert k30.weighted_geomean({"a": 0.25, "b": 1.0}, {"a": 1, "b": 1}) == pytest.approx(0.5, rel=1e-14)
    with pytest.raises(k30.Refused):
        k30.weighted_geomean({"a": 0.5}, {"a": 0.0})


@pytest.mark.parametrize("spoil,why", [
    (lambda s: s + "\n", "not at its pinned sha256"),
    (None, "lacks a column"),
    (lambda s: s.replace("0.7,0.75", "0.7,-0.75"), "not finite and positive"),
    (lambda s: s.replace("0.7,0.75", "0.7,0.0"), "not finite and positive"),
    (lambda s: s.replace("0.7,0.75", "0.7,nan"), "not finite and positive"),
    (lambda s: s.replace("0.7,0.75", "0.7,x"), "is not a number"),
    (lambda s: s.replace("Endo,grp,30", "Endo,grp,-30"), "wt_cells are not finite"),
    (lambda s: s.replace("Blood,grp", "Endo,grp"), "labels are empty or repeat"),
])
def test_the_table_is_refused_off_its_pin_or_not_well_formed(synthetic_table, monkeypatch, spoil, why):
    text_ = _table_text(extra_col=spoil is not None)
    text_ = text_.replace("rep2_depth_ratio", "rep3_depth_ratio") if spoil is None else spoil(text_)
    synthetic_table.write_text(text_)
    if "pinned" not in why:
        monkeypatch.setattr(k30, "TABLE_SHA256", _sha(synthetic_table))
    with pytest.raises(k30.Refused, match=why):
        k30.load_table()


def test_the_vectors_on_a_vocabulary_by_name():
    bs = harness._t3_head()
    t = k30.load_table()
    V = _vocab()
    vec = k30.depth_vectors(bs, V, t)
    for task in ("g1", "g2"):
        want = _want(task)
        for i, lab in enumerate(V.labels):  # by name, whatever the order; a label of no ratio (Extra) at 1.0
            r = want.get(lab, 1.0)
            assert vec[task]["dx1"][i] == r
            assert vec[task]["dxh"][i] == pytest.approx(r ** 0.5, rel=1e-15)
            assert vec[task]["dxg"][i] == t["global"][task]
    assert set(vec) == {"g1", "g2"} and all(set(v) == set(k30.DEPTH_ARMS) for v in vec.values())
    rec = k30.depth_record(V, t, vec, "g2")
    assert rec["pairing"] == {"g1": "rep2_depth_ratio", "g2": "rep1_depth_ratio"} and set(rec["vectors"]) == {"g2"}
    assert rec["table_labels_not_in_vocabulary"] == ["Gone"] and rec["listed_labels"] == {"g1": 8, "g2": 8}
    assert rec["vectors"]["g2"]["dx1"]["vCM1"] == 0.5 and rec["vectors"]["g2"]["dx1"]["Extra"] == 1.0
    assert k30.depth_why({"plain_g2/seed_0.json": json.loads(json.dumps(rec))}, ["plain_g2/seed_0.json"], t) == []


def test_a_listed_label_the_vocabulary_lacks_is_refused():
    bs = harness._t3_head()
    with pytest.raises(bs.Refused, match="pSHF"):  # pSHF is listed for both tasks
        k30.depth_vectors(bs, _vocab([x for x in LABELS if x != "pSHF"]), k30.load_table())
    k30.depth_vectors(bs, _vocab(list(LABELS)), k30.load_table())  # Gone (1.0 everywhere) is not needed


def test_the_distance_to_u3ds_other_replicate_ratio():
    bs = harness._t3_head()
    V = _vocab()
    vec = k30.depth_vectors(bs, V, k30.load_table())
    # each replicate's own ratios on V (rep2's column is what g1 takes, rep1's what g2 takes)
    own = {"rep1": np.array([_want("g2").get(t, 1.0) for t in V.labels]),
           "rep2": np.array([_want("g1").get(t, 1.0) for t in V.labels])}
    # U3d's S1 for the target rep1 (task g1) is ctx.targets[rep2].full: rep2's ratios, which g1's dx1 holds
    ctx = SimpleNamespace(targets=[SimpleNamespace(name=n, full=SimpleNamespace(ratio=own[n])) for n in own])
    assert np.array_equal(own["rep2"], vec["g1"]["dx1"]) and np.array_equal(own["rep1"], vec["g2"]["dx1"])
    assert k30.u3d_distance(ctx, vec) == {"g1": 0.0, "g2": 0.0}
    ctx.targets[1].full.ratio = own["rep2"] * 1.5
    assert k30.u3d_distance(ctx, vec) == {"g1": pytest.approx(0.5 * own["rep2"].max()), "g2": 0.0}


# ------------------------------------------------------------------------------------------------ the builds


def _pool(V, seed=5):
    bs = harness._t3_head()
    rng = np.random.default_rng(seed)
    C = rng.poisson(6.0, size=(N_POOL, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = np.arange(N_POOL) % V.L
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    return bs.t3d.Pool(C, lab, P, V.L), C, lab


@pytest.fixture
def fold(monkeypatch):
    """A synthetic reference whose vocabulary holds the scenario's states, K18's plain-world preparation stubbed as
    K18's and K29's tests stub it, the barred factor by hand; construct is the builder's own (t3_decompose.build)."""
    bs = harness._t3_head()
    V = _vocab()
    pool, C, lab = _pool(V)
    rules = {"components": {"stage": {"min_cells": 5}}}
    ref = bs.Reference(preset="gata4", rules=rules, genes=GENES, V=V, pool=pool,
                       coords=np.random.default_rng(1).normal(size=(N_POOL, 3)))
    rng = np.random.default_rng(7)
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, N_POOL)]
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    phi = np.exp(rng.normal(0, 0.5, G))
    file_genes = GENES[::-1] + ["extra"]
    barred_file = np.exp(rng.normal(0, 0.3, G + 1))
    phi_b = k18.in_gene_order(barred_file, file_genes, GENES)
    adata = ad.AnnData(X=bs.t3d.lognorm(C).astype(np.float32), var=pd.DataFrame(index=GENES))
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {}}
    monkeypatch.setattr(bs, "assemble", lambda r, built, panel, g: (built.X, None, {}))
    monkeypatch.setattr(k7, "_prepared", lambda ctx: prep)
    monkeypatch.setattr(harness, "count_for", lambda board, n: 60)
    monkeypatch.setattr(harness, "finish", lambda ctx, p, X, coords: X)
    monkeypatch.setattr(k18, "_BARRED", {"phi": barred_file, "genes": file_genes, "spec": SPEC, "record": {}})
    k30.load_depth()
    F = k30.fold_factors(bs, ref, per, phi, group, np.asarray(C.mean(0), dtype=float), phi_b, SPEC)

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": adata}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    return SimpleNamespace(bs=bs, ref=ref, V=V, F=F, ctx=ctx, phi_b=phi_b, C=C)


def _spy(monkeypatch):
    """Record, independently of the script's own observer, what t3_decompose.build hands rescale_depth and lognorm (the
    counts before the depth step and the cells as built)."""
    t3d = harness._t3_head().t3d
    seen = {"depth": [], "lognorm": []}
    real_rescale, real_lognorm = t3d.rescale_depth, t3d.lognorm

    def rescale(C, ratio, P, rng):
        seen["depth"].append({"C": np.asarray(C).copy(), "ratio": np.asarray(ratio).copy()})
        return real_rescale(C, ratio, P, rng)

    def lognorm(C):
        seen["lognorm"].append(np.asarray(C).copy())
        return real_lognorm(C)

    monkeypatch.setattr(t3d, "rescale_depth", rescale)
    monkeypatch.setattr(t3d, "lognorm", lognorm)
    return seen


def test_the_base_is_k18s_barred_form_value_for_value(fold):
    for seed in (0, 3):
        mine = k30.build_arm(fold.bs, fold.ref, k30.BASE, 60, seed, fold.F)
        theirs = k18.build_arm(fold.bs, fold.ref, k18.BARRED, 60, seed, fold.F)  # K18's own barred_form
        assert np.array_equal(mine.X, theirs.X) and np.array_equal(mine.row, theirs.row)


@pytest.mark.parametrize("task", ["g1", "g2"])
def test_each_depth_arm_differs_from_the_base_by_the_depth_step_alone(fold, monkeypatch, task):
    seen = _spy(monkeypatch)
    base = k30.build_arm(fold.bs, fold.ref, k30.BASE, 60, 4, fold.F)
    assert seen["depth"] == []  # the base has no depth step
    pre = seen["lognorm"][-1]  # the base's counts: rows drawn, factor applied
    for arm in k30.DEPTH_ARMS:
        with k30.task_set(task):
            built = k30.build_arm(fold.bs, fold.ref, arm, 60, 4, fold.F)
        step = seen["depth"][-1]
        assert np.array_equal(built.row, base.row) and np.array_equal(built.lab, base.lab)  # the same cells
        assert np.array_equal(step["C"], pre)  # the same factor draws: the counts before the step are the base's
        vec = fold.F["k30"]["ratio"][task][arm]
        assert np.array_equal(step["ratio"], vec[built.lab])  # each cell thinned by its label's ratio
        want = _want(task)
        assert np.allclose(vec, [want.get(t, 1.0) ** k30.EXPONENT[arm] if arm != "dxg" else
                                 k30.load_table()["global"][task] for t in fold.V.labels], rtol=1e-15, atol=0)
        assert not np.array_equal(built.X, base.X)
        assert seen["lognorm"][-1].sum() < pre.sum()  # most ratios are below 1: the arm is shallower than the base


def test_in_the_plain_world_every_arms_counts_before_the_step_are_barred_forms(fold, monkeypatch):
    """K18's plain build (its fold, its assemble stream) through this script's wrapper, seed by seed as the harness
    builds: g1's key gives rep2's ratios, g2's key rep1's, the base none; every depth arm's counts before the step are
    barred_form's built counts (the spy's and the script's own record agree), and the task is cleared after each
    build."""
    bs = fold.bs
    real, got = bs.construct, []

    def construct(r, n, g, **kw):
        got.append(None if kw.get("ratio") is None else np.asarray(kw["ratio"]).copy())
        return real(r, n, g, **kw)

    monkeypatch.setattr(bs, "construct", construct)
    seen = _spy(monkeypatch)
    with k30.k18_arms(k30.ARMS):
        for task in ("g1", "g2"):
            for seed in (1_000_002, 1_000_007):  # the harness's build seeds
                k30._K30["counts"] = {}
                seen["depth"].clear()
                seen["lognorm"].clear()
                for arm in k30.ARMS:
                    k30._build_plain(fold.ctx(task, seed), arm)
                    F = k18._FOLD["F"]
                    if arm == k30.BASE:
                        assert got[-1] is None
                    else:
                        assert np.array_equal(got[-1], F["k30"]["ratio"][task][arm])
                        if arm == "dx1":
                            assert got[-1].tolist() == [_want(task).get(t, 1.0) for t in fold.V.labels]
                    assert "task" not in k30._K30
                base_built = seen["lognorm"][0]
                assert len(seen["depth"]) == 3 and all(np.array_equal(s["C"], base_built) for s in seen["depth"])
                rec = k30.counts_record()
                assert rec["pre_step_equals_base"] == {"dx1": True, "dxh": True, "dxg": True}
                assert {f["seed"] for f in rec["arms"].values()} == {seed}
                b = rec["arms"][k30.BASE]["mean_counts"]
                assert b == pytest.approx(base_built.sum(1).mean())
                assert all(rec["arms"][a]["pre_step_mean_counts"] == b for a in k30.DEPTH_ARMS)
                assert all(rec["arms"][a]["mean_counts"] < b for a in k30.DEPTH_ARMS)
    assert k18.ARMS == k30._K18_ORIGINAL["ARMS"] and k18._FOLD == {}


def test_the_count_figure_hashes_dtype_shape_and_bytes_in_place():
    C = np.arange(12, dtype=np.int32).reshape(3, 4)
    f = k30._counts_figure(C)
    assert f == {"mean_counts": 22.0, "cells": 3, "sha256": k30._counts_figure(C.copy())["sha256"]}
    assert k30._counts_figure(C.astype(np.int64))["sha256"] != f["sha256"]  # another dtype
    assert k30._counts_figure(C.reshape(4, 3))["sha256"] != f["sha256"]  # another shape
    D = C.copy()
    D[2, 3] += 1
    assert k30._counts_figure(D)["sha256"] != f["sha256"]  # one count
    assert k30._counts_figure(np.asfortranarray(C))["sha256"] == f["sha256"]  # the layout does not matter
    assert k30._counts_figure(np.zeros((0, 4), dtype=np.int32)) == {
        "mean_counts": 0.0, "cells": 0, "sha256": hashlib.sha256(b"<i4(0, 4)").hexdigest()}


def test_the_count_record_flags_an_arm_whose_counts_before_the_step_differ():
    k30._K30.update(task="g1", counts={})
    k30._note_counts(k30.BASE, 3, {"built": {"mean_counts": 10.0, "cells": 5, "sha256": "a"}})
    k30._note_counts("dx1", 3, {"built": {"mean_counts": 6.0, "cells": 5, "sha256": "b"},
                                "pre_step": {"mean_counts": 10.0, "cells": 5, "sha256": "a"}})
    k30._note_counts("dxh", 3, {"built": {"mean_counts": 8.0, "cells": 5, "sha256": "c"},
                                "pre_step": {"mean_counts": 10.5, "cells": 5, "sha256": "z"}})
    k30._note_counts("dxg", 4, {"built": {"mean_counts": 7.0, "cells": 5, "sha256": "d"},
                                "pre_step": {"mean_counts": 10.0, "cells": 5, "sha256": "a"}})  # another seed
    assert k30.counts_record()["pre_step_equals_base"] == {"dx1": True, "dxh": False, "dxg": False}
    k30._K30.pop("task")
    k30._note_counts("dx1", 3, {"built": {"mean_counts": 1.0, "cells": 1, "sha256": "e"}})  # no task: not kept
    assert k30.counts_record()["arms"]["dx1"]["built_sha256"] == "b"


def test_a_depth_arm_needs_a_task_and_the_task_is_cleared(fold):
    with pytest.raises(RuntimeError, match="no task is set"):
        k30.build_arm(fold.bs, fold.ref, "dx1", 60, 1, fold.F)
    with pytest.raises(ValueError, match="boom"), k30.task_set("g2"):
        assert k30._K30["task"] == "g2"
        raise ValueError("boom")
    assert "task" not in k30._K30
    with pytest.raises(RuntimeError, match="unknown task"), k30.task_set("m"):
        pass
    with pytest.raises(RuntimeError, match="unknown arm"):
        k30.build_arm(fold.bs, fold.ref, "refused_form", 60, 1, fold.F)
    other = SimpleNamespace(**{**fold.ref.__dict__, "V": _vocab()})
    with pytest.raises(RuntimeError, match="another vocabulary"), k30.task_set("g1"):
        k30.build_arm(fold.bs, other, "dx1", 60, 1, fold.F)
    t3d = harness._t3_head().t3d
    assert t3d.rescale_depth.__name__ == "rescale_depth" and t3d.lognorm.__name__ == "lognorm"  # observers restored


def test_k18s_names_are_replaced_inside_the_context_and_restored():
    orig = (k18.ARMS, k18.build_arm, k18.fold_factors)
    with k30.k18_arms(k30.ARMS):
        assert k18.ARMS == k30.ARMS and k18.build_arm is k30.build_arm and k18.fold_factors is k30.fold_factors
    assert (k18.ARMS, k18.build_arm, k18.fold_factors) == orig
    with pytest.raises(ValueError, match="boom"), k30.k18_arms((k30.BASE,)):
        raise ValueError("boom")
    assert (k18.ARMS, k18.build_arm, k18.fold_factors) == orig
    for bad in (("dx1",), (k30.BASE, "own_d")):
        with pytest.raises(RuntimeError), k30.k18_arms(bad):
            pass
    assert k30.part_arms("guard") == ("barred_form",) and k30.part_arms("score") == k30.ARMS
    k30._K30.clear()
    with pytest.raises(RuntimeError, match="not loaded"):
        k30.fold_factors(None, None, {}, None, None, None, None, SPEC)


def test_in_the_dm_world_the_step_goes_on_k18s_matched_pool_by_u3ds_path(monkeypatch):
    """K18's dm_row on a synthetic context (the frozen draws and orth_diag's thinning for real; the estimate, scorer and
    de_genes stubbed, as K18's own test does), through this script's _dm_task: rep1's row is g1's and takes rep2's
    ratios, rep2's row is g2's and takes rep1's; the check arm and the base take none; every arm draws from the matched
    pool (U3d's path: t3_decompose.build with the matched pool as its source and the ratio as `ratio`), and every
    depth arm's counts before the step are barred_form's built counts."""
    from common import core_metrics

    from vec import evaluate

    bs = harness._t3_head()
    t3d = bs.t3d
    rng = np.random.default_rng(3)
    V = _vocab()
    pool, C, lab = _pool(V, seed=9)
    rules = {"protocol": {"working_copy_fraction": 0.5, "split_seed": 0, "prediction_cells_multiplier": 1},
             "components": {"stage": {"min_cells": 5}}}
    targets = [SimpleNamespace(name=n, X_rel=rng.normal(size=(40, G)).astype(np.float32), hash_rel=np.arange(40),
                               hash_pool=np.arange(60)) for n in ("rep1", "rep2")]
    ref_C = rng.poisson(20.0, size=(50, G)).astype(np.int32)
    ctx = SimpleNamespace(targets=targets, rules=rules, V=V, ref_X=t3d.lognorm(ref_C), preset_name="gata4",
                          genes=GENES, wt=SimpleNamespace(pool=pool), stages=[])
    group = np.array(["A", "B"])[rng.integers(0, 2, N_POOL)]
    per = {"A": np.exp(rng.normal(0, 0.3, G)), "B": np.exp(rng.normal(0, 0.3, G))}
    phi, phi_b = np.exp(rng.normal(0, 0.3, G)), np.exp(rng.normal(0, 0.3, G))
    refw = bs.Reference(preset="gata4", rules=rules, genes=GENES, V=V, pool=pool, coords=np.zeros((N_POOL, 3)))
    k30.load_depth()
    monkeypatch.setattr(t3d, "estimate", lambda c, t, rows, g: SimpleNamespace(ratio=np.full(V.L, 0.8)))
    monkeypatch.setattr(evaluate, "panel_metrics", lambda *a, **k: {})
    monkeypatch.setattr(core_metrics, "de_genes", lambda a, b: (np.arange(5), np.arange(5, 10), None))
    monkeypatch.setattr(k18, "_scored", lambda X, truth, ref, anchors, up, dn, seed: {
        "score": 50.0, "skills": {}, "severity": {"state": "scored", "beta": 0.1, "r2": 0.1}})
    seen = _spy(monkeypatch)
    real_build, calls = t3d.build, []

    def build(c, n, g, f, src=None, ratio=None, **kw):
        k = len(seen["depth"])
        out = real_build(c, n, g, f, src=src, ratio=ratio, **kw)
        calls.append({"ratio": None if ratio is None else np.asarray(ratio).copy(), "pool": c.wt.pool if src is None
                      else src, "pre": seen["depth"][k]["C"] if ratio is not None else None,
                      "built": seen["lognorm"][-1]})
        return out

    monkeypatch.setattr(t3d, "build", build)
    with k30.k18_arms(k30.ARMS):
        F = k18.fold_factors(bs, refw, per, phi, group, ctx.ref_X.mean(0), phi_b, SPEC)
        monkeypatch.setattr(k18, "_STATE", {"ctx": ctx, "phi_b": phi_b, "ref_counts": ref_C,
                                            "ref_lab": rng.integers(0, V.L, 50), "F": F, "threads": 1,
                                            "ref_props": ref_C / ref_C.sum(1, keepdims=True)})
        for ti, task in ((0, "g1"), (1, "g2")):
            calls.clear()
            row = k30._dm_task((ti, 4))
            assert row["target"] == ("rep1", "rep2")[ti] and set(row["scores"]) == set(k30.ARMS)
            check, base, *arms = calls
            assert check["ratio"] is None and base["ratio"] is None  # the check arm, then the base: no depth step
            matched = check["pool"]
            assert matched is not pool and matched.C.sum() < 0.9 * pool.C.sum()  # K18's matched pool, thinned
            for got, arm in zip(arms, k30.DEPTH_ARMS, strict=True):
                assert np.array_equal(got["ratio"], F["k30"]["ratio"][task][arm])
                assert got["pool"] is matched and base["pool"] is matched  # every arm on the matched pool
                assert np.array_equal(got["pre"], base["built"])  # the counts before the step are barred_form's
                assert got["built"].sum() < base["built"].sum()
            assert row["counts"]["pre_step_equals_base"] == {"dx1": True, "dxh": True, "dxg": True}
            assert "task" not in k30._K30
    rec = F["k30"]["record"]
    assert rec["vectors"]["g1"]["dx1"]["aCM1"] == 0.4 and rec["vectors"]["g2"]["dx1"]["vCM1"] == 0.5


# ------------------------------------------------------------------------------------------------ the guard


BASE_SCORE = 56.0


def _refs() -> dict:
    """K18's committed barred_form rows by hand."""
    return {tuple(c.split()): [{"seed": s, "scores": {k18.BARRED: BASE_SCORE + s}} for s in SEEDS] for c in k30.CELLS}


def _guard_rows(shift=None, one=None, drop=None) -> dict:
    out = {}
    for c in k30.CELLS:
        rows = [{"seed": s, "scores": {k30.BASE: BASE_SCORE + s + (shift or {}).get(c, 0.0)}} for s in SEEDS]
        if one and c in one:
            rows[3]["scores"][k30.BASE] += one[c]
        out[c] = rows[:9] if drop == c else rows
    return out


def test_the_void_guard_within_k17s_tolerance_in_all_four_cells():
    refs = _refs()
    ok = k30.guard_eval(_guard_rows(), refs)
    assert ok["passes"] and ok["failed"] == [] and set(ok["cells"]) == set(k30.CELLS)
    assert k30.guard_eval(_guard_rows(shift={"dm g2": 0.09}, one={"plain g1": 0.49}), refs)["passes"]
    assert k30.guard_eval(_guard_rows(shift={"plain g2": 0.0625}, one={"dm g1": -0.5 + 0.0625}), refs)["passes"]
    for kw in ({"shift": {"plain g2": 0.11}}, {"shift": {"dm g1": -0.11}}, {"one": {"dm g2": 0.51}},
               {"one": {"plain g1": -0.6}}, {"drop": "plain g1"}):
        res = k30.guard_eval(_guard_rows(**kw), refs)
        assert not res["passes"] and res["failed"] == ["barred_form"], kw
        assert sum(not v["passes"] for v in res["cells"].values()) == 1
    res = k30.guard_eval(_guard_rows(), {k: v for k, v in refs.items() if k != ("dm", "g1")})
    assert not res["passes"] and not res["cells"]["dm g1"]["passes"]  # every cell is needed


def test_the_run_steps_base_is_matched_with_the_guard_steps():
    g = _guard_rows()
    assert k30.match_eval(_guard_rows(shift={"dm g1": 0.09}), g)["passes"]
    for kw in ({"shift": {"plain g1": 0.11}}, {"one": {"dm g2": -0.6}}, {"drop": "dm g1"}):
        res = k30.match_eval(_guard_rows(**kw), g)
        assert not res["passes"] and sum(not v["passes"] for v in res["cells"].values()) == 1, kw


# ------------------------------------------------------------------------------------------------ the rule


def _state(gated, arm: str, seed: int) -> str:
    g = (gated or {}).get(arm, ())
    if isinstance(g, dict):
        return g.get(seed, "scored")
    return "r2_gate" if seed in g else "scored"


MEAN_COUNTS = {k30.BASE: 128.0, "dx1": 72.0, "dxh": 96.0, "dxg": 64.0}


def _counts(arms) -> dict:
    """A count record by hand: each arm's mean counts per built cell, every depth arm's counts before the step the
    base's."""
    return {"arms": {a: {"seed": 0, "mean_counts": MEAN_COUNTS[a]} for a in arms if a in MEAN_COUNTS},
            "pre_step_equals_base": {a: True for a in arms if a in k30.DEPTH_ARMS}}


def _row_set(table: dict[str, list[float]], gated=None) -> list[dict]:
    out = []
    for i, s in enumerate(SEEDS):
        sc = {a: float(v[i]) for a, v in table.items()}
        out.append({"seed": s, "scores": sc, "counts": _counts(sc),
                    "skills": {a: {"de_score": sc[a] / 128.0, "mmd_u": 0.5} for a in sc},
                    "t3": {a: {"severity": {"state": _state(gated, a, s), "beta": 0.1, "r2": 0.05}} for a in sc}})
    return out


def _cell_table(c: str, offsets: dict, base=lambda s: BASE_SCORE + s) -> dict:
    """The base at `base(seed)`; each depth arm the base plus its offset (a number, or ten numbers by seed); the floor
    in the plain cells. Offsets are dyadic or chosen so that the paired means are exact."""
    t = {k30.BASE: [base(s) for s in SEEDS]}
    for a in k30.DEPTH_ARMS:
        v = offsets.get(a, 0.0)
        t[a] = [base(s) + x for s, x in zip(SEEDS, v if isinstance(v, list) else [v] * len(SEEDS), strict=True)]
    if c.startswith("plain"):
        t[harness.FLOOR] = [50.0 for _ in SEEDS]
    return t


GOOD = {"plain g1": {"dx1": 8.0, "dxh": 4.0, "dxg": 6.0}, "plain g2": {"dx1": 8.0, "dxh": 4.0, "dxg": 6.0},
        "dm g1": {"dx1": -4.0, "dxh": -1.0, "dxg": -5.0}, "dm g2": {"dx1": -4.0, "dxh": -1.0, "dxg": -5.0}}
# EV at plain: dx1 8, dxh 4, dxg 6; at 0.5: dx1 2.0, dxh 1.5, dxg 0.5: dx1 is named at both


def _per(cells=None, gated=None, base=lambda s: BASE_SCORE + s) -> dict:
    per = {}
    for c in k30.CELLS:
        t = _cell_table(c, (cells or GOOD).get(c, {}), base)
        per[c] = k30.summarise(_row_set(t, (gated or {}).get(c)), tuple(t))
    return per


def test_ev_is_exact_at_each_weight():
    assert k30.ev([8.0, 8.0], [-4.0, -4.0], "plain") == 8 and k30.ev([8.0, 8.0], [-4.0, -4.0], "0.5") == 2
    assert k30.ev([1.75, 1.75], [0.25, 0.25], "0.5") == 1  # exactly 1
    assert k30.ev([2.0, 0.0], [1.0, -1.0], "0.5") == Fraction(1, 2)  # the means of the two plain and two dm cells
    assert k30.ev([0.3, 0.3], [-9.0, 7.0], "plain") == Fraction(0.3)  # the dm cells weigh nothing at plain
    with pytest.raises(KeyError):
        k30.ev([1.0, 1.0], [1.0, 1.0], "0.7")
    with pytest.raises(ValueError, match="unknown S3 weight"):
        k30.decide(_per(), "0.7")


@pytest.mark.parametrize("weight,evs", [("plain", {"dx1": 8.0, "dxh": 4.0, "dxg": 6.0}),
                                        ("0.5", {"dx1": 2.0, "dxh": 1.5, "dxg": 0.5})])
def test_a_named_arm_and_the_table(weight, evs):
    dec = k30.decide(_per(), weight)
    assert dec["clauses"] == [] and dec["named"] == "dx1" and dec["file"] == "T3_c2_s2_dd.h5ad"
    assert dec["outcome"] == ("Try K30: NAME dx1: the file may be built twice and uploaded under the captain's bound "
                              "yes (R53).")
    assert {a: dec["arms"][a]["ev"] for a in k30.DEPTH_ARMS} == evs and dec["s3_weight"] == weight
    assert dec["arms"]["dx1"]["plain_mean"] == 8.0 and dec["arms"]["dx1"]["dm_mean"] == -4.0
    assert all(dec["arms"][a]["may_be_named"] for a in k30.DEPTH_ARMS)
    assert all(dec["table"]["dx1"][c]["AHEAD"] for c in PLAIN) and not any(dec["table"]["dx1"][c]["AHEAD"] for c in DM)
    assert dec["stop"] == {"S1": False, "S2": False, "S3": False, "S4": False} and dec["ev_tied_with_best"] == []
    lines = k30.table_lines(dec)
    assert any(x.startswith("  dx1: plain g1 +8.00 (sd 0.00, 10 of 10 positive) AHEAD") for x in lines)
    assert lines[-1].startswith("  severity-gated seeds of 10 per arm and cell: plain g1 barred_form 0, dx1 0")
    assert all(c in lines[-1] for c in k30.CELLS)


def test_the_weight_can_change_the_name():
    cells = {c: {"dx1": 8.0, "dxh": 4.0} for c in PLAIN} | {c: {"dx1": -10.0, "dxh": -1.0} for c in DM}
    assert k30.decide(_per(cells), "plain")["named"] == "dx1"  # EV 8 against 4
    dec = k30.decide(_per(cells), "0.5")  # EV -1.0 against 1.5
    assert dec["named"] == "dxh" and dec["arms"]["dx1"]["ev"] == -1.0


def test_at_plain_the_dm_cells_decide_nothing():
    plain = {c: {"dx1": 3.0, "dxh": 2.0, "dxg": 1.5} for c in PLAIN}
    outs = set()
    for dm in ({}, {"dx1": -50.0}, {"dx1": 40.0, "dxh": -40.0, "dxg": 9.0}, {"dxh": 60.0}):
        for gated in (None, {"dm g1": {"dx1": list(range(10))}, "dm g2": {"dxh": list(range(10))}}):
            dec = k30.decide(_per(plain | {c: dm for c in DM}, gated=gated), "plain")
            outs.add((dec["outcome"], tuple(dec["clauses"]), tuple((a, dec["arms"][a]["ev"]) for a in k30.DEPTH_ARMS)))
    assert outs == {(k30.NAME_LINE.format(arm="dx1"), (), (("dx1", 3.0), ("dxh", 2.0), ("dxg", 1.5)))}


def test_the_case_where_nothing_changes_is_stop_s1_and_s4():
    for weight in ("plain", "0.5"):
        dec = k30.decide(_per({c: {} for c in k30.CELLS}), weight)
        assert dec["clauses"] == ["S1", "S4"] and dec["named"] is None
        assert dec["outcome"] == "Try K30: STOP (S1, S4): no file; C2 stays."
        assert all(dec["arms"][a]["ev"] == 0.0 for a in k30.DEPTH_ARMS)


@pytest.mark.parametrize("plain,clauses,named", [
    ({"dx1": 1.0}, ["S1"], None),  # a mean of exactly +1.0 is not above +1.0
    ({"dx1": [2.0] * 7 + [-0.25] * 3}, ["S1"], None),  # mean 1.325, 7 of 10 positive
    ({"dx1": [2.0] * 7 + [0.0] * 3}, ["S1"], None),  # a seed at exactly 0 is not positive
    ({"dx1": [2.0] * 8 + [-0.5] * 2}, [], "dx1"),  # mean 1.5 with exactly 8 of 10 positive: AHEAD
])
def test_ahead_at_its_edges(plain, clauses, named):
    cells = {c: plain for c in PLAIN} | {c: {} for c in DM}
    dec = k30.decide(_per(cells), "plain")
    assert dec["clauses"] == clauses and dec["named"] == named


def test_s1_needs_one_arm_ahead_in_both_plain_cells():
    cells = {"plain g1": {"dx1": 3.0, "dxh": 0.0}, "plain g2": {"dx1": 0.0, "dxh": 3.0}, "dm g1": {}, "dm g2": {}}
    dec = k30.decide(_per(cells), "plain")  # each AHEAD in one plain cell only: no arm passes S1
    assert dec["clauses"] == ["S1"] and not any(dec["arms"][a]["passes_s1"] for a in k30.DEPTH_ARMS)
    dm_only = {"plain g1": {}, "plain g2": {}, "dm g1": {"dx1": 9.0}, "dm g2": {"dx1": 9.0}}
    assert k30.decide(_per(dm_only), "0.5")["clauses"] == ["S1", "S4"]  # AHEAD in the dm cells does not pass S1


@pytest.mark.parametrize("gated,s2,named", [
    ({"plain g1": {"dx1": [0, 1]}}, True, None),  # 2 > 0 + 1 in a plain cell, the only arm that passes S1
    ({"plain g2": {"dx1": [0, 1, 2], "barred_form": [5]}}, True, None),  # 3 > 1 + 1
    ({"plain g1": {"dx1": [0]}}, False, "dx1"),  # 1 is not more than 0 + 1
    ({"plain g2": {"dx1": [0, 1], "barred_form": [7]}}, False, "dx1"),  # 2 is not more than 1 + 1
    ({"plain g1": {"dx1": {0: "beta_gate", 1: "no_response"}}}, True, None),  # every gated state counts
    ({"dm g1": {"dx1": [0, 1, 2, 3]}}, False, "dx1"),  # S2 reads the plain cells only
    ({"plain g1": {"barred_form": [0, 1, 2]}}, False, "dx1"),  # fewer gated than the base never stops
])
def test_s2_counts_gated_seeds_in_the_plain_cells_against_the_base_plus_one(gated, s2, named):
    only_dx1 = {c: {"dx1": 3.0} for c in PLAIN} | {c: {} for c in DM}  # dx1 alone passes S1
    dec = k30.decide(_per(only_dx1, gated=gated), "plain")
    assert dec["stop"]["S2"] is s2 and dec["named"] == named
    assert dec["clauses"] == (["S2"] if s2 else []) and bool(dec["arms"]["dx1"]["s2_cells"]) is s2
    for c in k30.CELLS:  # gated seeds are reported per arm and cell, all four cells
        for a in k30.ARMS:
            g = (gated.get(c) or {}).get(a, ())
            assert dec["gated"][c][a] == len(g)


def test_s2_falls_back_to_the_next_arm_that_passes_s1_and_s2():
    cells = {c: {"dx1": 3.0, "dxh": 2.0} for c in PLAIN} | {c: {} for c in DM}
    dec = k30.decide(_per(cells, gated={"plain g2": {"dx1": [0, 1, 2]}}), "plain")
    assert dec["clauses"] == [] and dec["named"] == "dxh" and not dec["arms"]["dx1"]["may_be_named"]
    assert dec["arms"]["dx1"]["s2_cells"] == ["plain g2"] and dec["best"] == "dxh"


def test_s3_reads_the_best_ev_among_the_arms_that_may_be_named():
    exact = {c: {"dx1": 1.75} for c in PLAIN} | {c: {"dx1": 0.25} for c in DM}  # EV at 0.5 exactly +1.0
    dec = k30.decide(_per(exact), "0.5")
    assert dec["clauses"] == ["S3"] and dec["named"] is None and dec["arms"]["dx1"]["ev_exact"] == "1/1"
    assert k30.decide(_per(exact), "plain")["named"] == "dx1"  # at plain the EV is the plain mean, 1.75
    above = {c: {"dx1": 1.75} for c in PLAIN} | {c: {"dx1": 0.5} for c in DM}  # EV 1.125
    assert k30.decide(_per(above), "0.5")["named"] == "dx1"
    # dx1's EV 1.0 at 0.5 is not to be named (S2); dxh's EV 0.25 is the best that may be: S3
    cells = {c: {"dx1": 3.0, "dxh": 1.25} for c in PLAIN} | {c: {"dx1": -1.0, "dxh": -0.75} for c in DM}
    dec = k30.decide(_per(cells, gated={"plain g1": {"dx1": [0, 1, 2]}}), "0.5")
    assert dec["clauses"] == ["S3"] and dec["best"] == "dxh" and dec["arms"]["dxh"]["ev"] == 0.25


def test_at_plain_s3_never_stops_an_arm_that_passes_s1():
    for p in (1.0625, 1.5, 2.0):
        cells = {c: {"dxg": p} for c in PLAIN} | {c: {"dxg": -100.0} for c in DM}
        dec = k30.decide(_per(cells), "plain")
        assert not dec["stop"]["S3"] and dec["named"] == "dxg"


def test_s4_at_its_edge():
    zero = lambda s: 0.0  # noqa: E731  (a base of 0: each paired difference is the arm's score exactly)
    edge = {"plain g1": {"dxg": [-3.0] + [0.0] * 9}, "plain g2": {}, "dm g1": {}, "dm g2": {}}  # |mean| exactly 0.3
    dec = k30.decide(_per(edge, base=zero), "plain")
    assert dec["table"]["dxg"]["plain g1"]["mean"] == -0.3 and dec["clauses"] == ["S1"]
    inside = {"plain g1": {"dxg": [2.875] + [0.0] * 9}, "plain g2": {"dx1": -0.25}, "dm g1": {"dx1": -5.0}}
    assert k30.decide(_per(inside | {"dm g2": {}}, base=zero), "plain")["clauses"] == ["S1", "S4"]


def test_ties_in_ev_go_to_the_first_of_dx1_dxh_dxg():
    # at plain, dx1 and dxh tie on their plain means; the dm cells differ and decide nothing: dx1
    cells = {c: {"dx1": 2.0, "dxh": 2.0} for c in PLAIN} | {c: {"dx1": -1.0, "dxh": 5.0} for c in DM}
    dec = k30.decide(_per(cells), "plain")
    assert dec["named"] == "dx1" and dec["ev_tied_with_best"] == ["dxh"]
    # at 0.5: dx1 plain 4.0, dm -1.0, EV 3/2; dxh plain 2.25, dm 0.75, EV 3/2 exactly; both may be named: dx1
    cells = {c: {"dx1": 4.0, "dxh": 2.25} for c in PLAIN} | {c: {"dx1": -1.0, "dxh": 0.75} for c in DM}
    dec = k30.decide(_per(cells), "0.5")
    assert dec["arms"]["dx1"]["ev_exact"] == dec["arms"]["dxh"]["ev_exact"] == "3/2"
    assert dec["named"] == "dx1" and dec["ev_tied_with_best"] == ["dxh"]
    same = {c: {"dxh": 2.0, "dxg": 2.0} for c in PLAIN} | {c: {"dxh": -1.0, "dxg": -1.0} for c in DM}
    assert k30.decide(_per(same), "plain")["named"] == "dxh"


def test_ev_mixes_the_worlds_at_one_half():
    cells = {c: {"dx1": 8.0, "dxg": 6.0} for c in PLAIN} | {c: {"dx1": -8.0, "dxg": -1.0} for c in DM}
    dec = k30.decide(_per(cells), "0.5")
    assert dec["arms"]["dx1"]["ev"] == 0.0 and dec["arms"]["dxg"]["ev"] == 2.5 and dec["named"] == "dxg"
    assert k30.decide(_per(cells), "plain")["named"] == "dx1"


def test_summarise_reports_without_deciding():
    per = _per()
    assert set(per["plain g1"]["paired"]) == {"dx1 - barred_form", "dxh - barred_form", "dxg - barred_form",
                                              "dx1 - dxh", "dx1 - dxg", "barred_form - floor"}
    assert "barred_form - floor" not in per["dm g1"]["paired"]
    assert per["plain g2"]["paired"]["dx1 - dxh"]["mean"] == 4.0
    assert per["plain g2"]["per_term"]["dx1 - barred_form"]["de_score"] == pytest.approx(8.0 / 128)
    assert per["dm g1"]["counts"]["dx1"] == {"mean_counts": 72.0, "against_base": 72.0 / 128.0,
                                             "pre_step_equals_base_seeds": 10}
    assert per["dm g1"]["counts"]["barred_form"]["pre_step_equals_base_seeds"] is None
    assert "floor" not in per["plain g1"]["counts"]
    per = _per(gated={"dm g2": {"dxh": [3, 4]}})
    assert per["dm g2"]["flips"] == {"dx1": [], "dxh": [3, 4], "dxg": []}
    assert per["dm g2"]["arms"]["dxh"]["severity_gated"] == 2


# ------------------------------------------------------------------------------------------------ the reading


def _prov(commit="abc", clean=True, machine="m1") -> dict:
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "x"}


@pytest.fixture
def records():
    """Each task's depth record as the runs write it, on a synthetic vocabulary, from the synthetic table."""
    bs = harness._t3_head()
    t = k30.load_table()
    V = _vocab()
    vec = k30.depth_vectors(bs, V, t)
    return {n: json.loads(json.dumps(k30.depth_record(V, t, vec, n))) for n in k30.TASKS}


def _body(world: str, task: str, row: dict, part: str, arms, prov: dict, rec: dict, weight: str) -> dict:
    counts = row.pop("counts")
    if world == "plain":
        return {"task": {"key": k30.TASKS[task].key}, "rows": [row], "seed": row["seed"], "world": "plain",
                "part": part, "arms": list(arms), "depth": rec, "counts": counts, "s3_weight": weight,
                "provenance": prov}
    return {**row, "target": k18.DM_TARGETS[task], "task": task, "world": "dm", "part": part, "arms": list(arms),
            "depth": rec, "counts": counts, "u3d_other_ratio_largest_distance": 0.0, "s3_weight": weight,
            "provenance": prov}


def _write_run(root: Path, records: dict, cells=None, gcells=None, commit="abc", clean=True, machine="m1",
               guard_passes=True, score=True, bad=None, gated=None, weight="plain") -> SimpleNamespace:
    """A run folder by hand (the guard part, guard.json, run_started.json and the scoring part) and K18's committed
    files by hand beside it, their pins computed here. `cells` and `gcells` give the offsets per cell, of each depth
    arm (and of the base as 'barred_form', against K18's rows); `bad` maps a file's path to a function that spoils its
    body (None: the file is not written)."""
    cells = GOOD if cells is None else cells
    root.mkdir(parents=True)
    prov = _prov(commit, clean, machine)
    gsha = {}
    for part, prefix, arms in (("guard", "guard/", (k30.BASE,)), ("score", "", k30.ARMS)):
        if part == "score" and not score:
            continue
        src = (cells if gcells is None else gcells) if part == "guard" else cells
        for c in k30.CELLS:
            world, task = c.split()
            off = src.get(c, {})
            shift = off.get(k30.BASE, 0.0)
            t = _cell_table(c, off, base=lambda s, shift=shift: BASE_SCORE + s + (
                shift[s] if isinstance(shift, list) else shift))
            want = ((harness.FLOOR,) if world == "plain" else ()) + tuple(arms)
            for r in _row_set({a: t[a] for a in want}, (gated or {}).get(c) if part == "score" else None):
                rel = prefix + (f"plain_{task}/seed_{r['seed']}.json" if world == "plain" else
                                f"dm/{task}_seed_{r['seed']}.json")
                body = _body(world, task, r, part, arms, prov, records[task], weight)
                if bad and rel in bad:
                    body = bad[rel](body)
                if body is None:
                    continue
                p = root / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                k18.write_task_file(p, body)
                if part == "guard":
                    gsha[rel] = _sha(p)
        if part == "guard":
            k18.write_task_file(root / "guard.json", {"passes": guard_passes,
                                                      "failed": [] if guard_passes else ["barred_form"], "cells": {},
                                                      "files_sha256": gsha, "provenance": prov})
        else:
            (root / "run_started.json").write_text(json.dumps({"row": "K30", "s3_weight": weight, "commit": commit,
                                                               "utc": "x"}))
    refs = root.parent / f"k18_{root.name}"
    refs.mkdir()
    rows = _refs()
    for n in k30.TASKS:
        (refs / f"plain_{n}.json").write_text(json.dumps({"rows": rows[("plain", n)]}))
    (refs / "dm.json").write_text(json.dumps({"rows": {n: rows[("dm", n)] for n in k30.TASKS}}))
    return SimpleNamespace(root=root, refs=refs, pins={f: _sha(refs / f) for f in k30.K18_SHA256}, k18rows=rows)


def _read(w, monkeypatch, weight="plain") -> dict:
    monkeypatch.setattr(k30, "K18_SHA256", w.pins)
    return k30.reading(w.root, w.refs, weight=weight)


def _drop(arm: str):
    def spoil(body):
        row = body["rows"][0] if "rows" in body else body
        for k in ("scores", "skills", "t3"):
            row[k].pop(arm)
        return body

    return spoil


def _set(key: str, **kw):
    def spoil(body):
        body[key] = {**body[key], **kw}
        return body

    return spoil


def test_the_reading_of_a_run_folder_that_names_an_arm(tmp_path, monkeypatch, records):
    res = _read(_write_run(tmp_path / "ok", records), monkeypatch)
    assert not res["void"], res.get("why")
    assert res["outcome"] == ("Try K30: NAME dx1: the file may be built twice and uploaded under the captain's bound "
                              "yes (R53).")
    assert res["tries_counted"] == 3 and res["counted_arms"] == list(k30.COUNTED) and res["run_step_started"]
    assert res["guard"]["passes"] and res["guard_match"]["passes"] and res["s3_weight"] == "plain"
    assert res["commit"] == "abc" and res["machine"] == "m1" and res["run_started"]["row"] == "K30"
    assert len(res["files_sha256"]) == 1 + 40 + 40 and res["depth"]["pairing"] == k30.COLUMN
    assert res["u3d_other_ratio_largest_distance"] == {"dm g1": 0.0, "dm g2": 0.0}
    txt = k30.text(res)
    assert txt.splitlines()[0] == res["outcome"]
    for x in ("stop: S1 no, S2 no, S3 no, S4 no", "named: dx1 (R53 names its file T3_c2_s2_dd.h5ad), EV +8.000 at "
              "the S3 weight plain", "the plain cells' mean +8.00", "the dm cells' mean -4.00 (at its reference's "
              "depth) with severity gated in 0 and 0 of 10 seeds (barred_form 0 and 0)", "dx1 - dxg",
              "barred_form - floor", "per-term skill, dxh - barred_form", "gate flips against the base",
              "tries counted: 3 (dx1, dxh, dxg); the run step started: True; S3 weight plain; ruling row R53",
              "pairing g1 takes rep2_depth_ratio, g2 takes rep1_depth_ratio", "EV at the S3 weight plain = "
              "mean(plain), the dm cells deciding nothing", "severity-gated seeds of 10 per arm and cell",
              "dx1 72.0 (0.562; 10 of 10)", "barred_form 128.0 (1.000)", "U3d's other-replicate ratio dm g1 0.000"):
        assert x in txt, x


def test_a_stop_reading(tmp_path, monkeypatch, records):
    res = _read(_write_run(tmp_path / "s", records, cells={c: {} for c in k30.CELLS}), monkeypatch)
    assert not res["void"] and res["outcome"] == "Try K30: STOP (S1, S4): no file; C2 stays."
    assert res["tries_counted"] == 3 and k30.text(res).splitlines()[0] == res["outcome"]
    assert "S4: nothing changes" in k30.text(res)


@pytest.mark.parametrize("kw,why,tries", [
    ({"bad": {"plain_g1/seed_3.json": _drop("dxh")}}, "plain_g1/seed_3.json: missing arms ['dxh']", 3),
    ({"bad": {"dm/g2_seed_7.json": _drop("barred_form")}}, "dm/g2_seed_7.json: missing arms ['barred_form']", 3),
    ({"bad": {"plain_g2/seed_0.json": _drop("floor")}}, "plain_g2/seed_0.json: missing arms ['floor']", 3),
    ({"bad": {"guard/dm/g1_seed_5.json": _drop("barred_form")}}, "guard/dm/g1_seed_5.json: missing arms", 3),
    ({"bad": {"dm/g1_seed_4.json": _set("provenance", commit="def")}}, "results from different commits", 3),
    ({"bad": {"plain_g1/seed_9.json": _set("provenance", machine="m2")}}, "results from different machines", 3),
    ({"bad": {"dm/g2_seed_1.json": _set("provenance", tree_clean=False)}}, "dm/g2_seed_1.json: a dirty tree", 3),
    ({"clean": False}, "guard.json: a dirty tree", 3),
    ({"bad": {"dm/g2_seed_2.json": _set("depth", table_sha256="0" * 64)}}, "its depth record is not the pinned", 3),
    ({"bad": {"plain_g2/seed_6.json": _set("depth", pairing={"g1": "rep1_depth_ratio", "g2": "rep2_depth_ratio"})}},
     "plain_g2/seed_6.json: its depth record is not the pinned table's", 3),
    ({"bad": {"guard/plain_g1/seed_2.json": _set("depth", **{"global": {"g1": 0.5, "g2": 0.5}})}},
     "guard/plain_g1/seed_2.json: its depth record is not the pinned table's", 3),
    ({"bad": {"dm/g1_seed_8.json": lambda b: {**b, "depth": None}}}, "dm/g1_seed_8.json: its depth record", 3),
    ({"bad": {"plain_g2/seed_4.json": lambda b: {**b, "s3_weight": "0.5"}}},
     "1 task file(s) were run under another S3 weight than plain, e.g. plain_g2/seed_4.json", 3),
    ({"weight": "0.5"}, "80 task file(s) were run under another S3 weight than plain", 3),
    ({"gcells": {"dm g2": {k30.BASE: 0.2}}}, "the guard step's barred_form, read again, is off", 4),
    ({"cells": {**GOOD, "plain g2": {**GOOD["plain g2"], k30.BASE: 0.11}}, "gcells": {}},
     "the run step's barred_form is off the guard step's in plain g2", 4),
    ({"cells": {**GOOD, "dm g1": {**GOOD["dm g1"], k30.BASE: [0.0] * 3 + [0.6] + [0.0] * 6}}, "gcells": {}},
     "the run step's barred_form is off the guard step's in dm g1", 4),
])
def test_each_void_condition(tmp_path, monkeypatch, records, kw, why, tries):
    res = _read(_write_run(tmp_path / "v", records, **kw), monkeypatch)
    assert res["void"] and any(why in x for x in res["why"]), res["why"]
    assert res["outcome"].startswith("Try K30: VOID: ") and res["outcome"].endswith("; no file.")
    assert why in res["outcome"] and "decision" not in res
    assert res["tries_counted"] == tries  # three; four when the base fails its guard in the run step
    assert res["counted_arms"] == [*k30.COUNTED, *(["barred_form"] if tries == 4 else [])]
    assert "VOID (checked first; no STOP or NAME clause is applied)" in k30.text(res)


def test_a_depth_vector_off_the_pinned_table_voids(tmp_path, monkeypatch, records):
    def spoil(body):
        body["depth"] = json.loads(json.dumps(body["depth"]))  # a copy: the other files keep their record
        body["depth"]["vectors"]["g1"]["dxh"]["vCM1"] = 0.5  # vCM1 is at 1.0 in rep2's column (g1): its root is 1.0
        return body

    res = _read(_write_run(tmp_path / "d", records, bad={"plain_g1/seed_1.json": spoil}), monkeypatch)
    assert res["void"] and res["why"] == ["plain_g1/seed_1.json: its depth vectors are not the other replicate's "
                                          "ratios of the pinned table"]

    def swap(body):  # g2's file carrying g1's vectors: the pairing reversed
        body["depth"] = {**body["depth"], "vectors": {"g2": records["g1"]["vectors"]["g1"]}}
        return body

    res = _read(_write_run(tmp_path / "e", records, bad={"dm/g2_seed_0.json": swap}), monkeypatch)
    assert res["void"] and res["why"] == ["dm/g2_seed_0.json: its depth vectors are not the other replicate's ratios "
                                          "of the pinned table"]


def test_the_table_off_its_pin_at_read_time_voids(tmp_path, monkeypatch, records, synthetic_table):
    w = _write_run(tmp_path / "t", records)
    synthetic_table.write_text(_table_text() + "\n")
    res = _read(w, monkeypatch)
    assert res["void"] and "the depth table is not at its pinned sha256" in res["outcome"]


def test_k18s_committed_files_off_their_pins_void(tmp_path, monkeypatch, records):
    w = _write_run(tmp_path / "p", records)
    (w.refs / "dm.json").write_text((w.refs / "dm.json").read_text() + " ")
    res = _read(w, monkeypatch)
    assert res["void"] and "K18's committed dm.json is missing or not at its pinned sha256" in res["outcome"]
    assert res["tries_counted"] == 3 and res["guard"] is None


def test_a_guard_that_failed_in_the_guard_step_is_void_with_one_counted_try(tmp_path, monkeypatch, records):
    res = _read(_write_run(tmp_path / "g", records, guard_passes=False, score=False), monkeypatch)
    assert res["void"] and res["tries_counted"] == 1 and res["counted_arms"] == ["barred_form"]
    assert not res["run_step_started"]
    assert res["outcome"] == ("Try K30: VOID: the guard failed in the guard step: barred_form is off K18's committed "
                              "barred_form; no file.")
    assert k30.text(res).splitlines()[0] == res["outcome"]


@pytest.mark.parametrize("kw", [{"clean": False}, {"bad": {"dm/g1_seed_2.json": _drop("dx1")}},
                                {"bad": {"plain_g1/seed_1.json": _set("provenance", commit="later")}},
                                {"guard_passes": False}])
def test_void_is_checked_first_and_wins(tmp_path, monkeypatch, records, kw):
    assert _read(_write_run(tmp_path / "n", records), monkeypatch)["outcome"].startswith("Try K30: NAME dx1")
    res = _read(_write_run(tmp_path / "v", records, **kw), monkeypatch)
    assert res["void"] and res["outcome"].startswith("Try K30: VOID: ") and "decision" not in res
    assert "NAME" not in res["outcome"] and "STOP" not in res["outcome"]


def test_missing_files_are_not_a_reading_and_no_seed_file_is_opened(tmp_path, monkeypatch, records):
    w = _write_run(tmp_path / "i", records, bad={"plain_g2/seed_7.json": lambda b: None})
    (w.root / "dm" / "g1_seed_0.json").write_text("{ not JSON: the reading must not open this file")
    assert _read(w, monkeypatch) == {"incomplete": True, "void": False, "missing": ["plain_g2/seed_7.json"],
                                     "outcome": k30.INCOMPLETE}
    w2 = _write_run(tmp_path / "j", records, bad={"guard/dm/g1_seed_2.json": lambda b: None})
    (w2.root / "guard.json").write_text("{ not JSON: not opened while a guard file is missing")
    assert _read(w2, monkeypatch)["missing"] == ["guard/dm/g1_seed_2.json"]
    (w2.root / "guard.json").unlink()
    assert _read(w2, monkeypatch)["missing"] == ["guard.json", "guard/dm/g1_seed_2.json"]
    assert k30.text(_read(w2, monkeypatch)).startswith("Try K30: INCOMPLETE: not a reading")
    assert "At or after 2026-11-25 13:00 UTC a read writes the VOID reading 'INCOMPLETE past its cutoff'" in \
        k30.INCOMPLETE


def test_past_the_cutoff_a_missing_file_makes_the_void_reading(tmp_path, monkeypatch, records):
    monkeypatch.setattr(k30, "_now", lambda: k30.CUTOFF)  # at the cutoff itself
    w = _write_run(tmp_path / "late", records, bad={"plain_g2/seed_7.json": lambda b: None})
    res = _read(w, monkeypatch)
    assert res["void"] and res["late"] and res["tries_counted"] == 3 and res["run_step_started"]
    assert res["outcome"] == ("Try K30: VOID: INCOMPLETE past its cutoff (2026-11-25 13:00 UTC): 1 file(s) missing, "
                              "e.g. plain_g2/seed_7.json; no file.")
    w2 = _write_run(tmp_path / "early", records, score=False)  # the guard step done, the run step never started
    (w2.root / "guard.json").unlink()
    res = _read(w2, monkeypatch)
    assert res["void"] and res["tries_counted"] == 0 and res["outcome"].startswith("Try K30: VOID: INCOMPLETE past")
    monkeypatch.setattr(k30, "_now", lambda: datetime(2026, 11, 25, 12, 59, 59, tzinfo=UTC))
    assert _read(w, monkeypatch)["incomplete"]


# ------------------------------------------------------------------------------------------------ the steps


@pytest.fixture
def ready(tmp_path, monkeypatch, records):
    """Rows K30 and R53 in a hand-made ledger, a clean stamp at commit abc on m1, a hand-made run folder and K18's
    files."""
    tries = tmp_path / "TRIES.md"
    tries.write_text(LEDGER)
    monkeypatch.setattr(k30, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: _prov())
    w = _write_run(tmp_path / "run", records)
    monkeypatch.setattr(k30, "K18_SHA256", w.pins)
    monkeypatch.setattr(k30, "reading", lambda runs, f=k30.reading, **kw: f(runs, w.refs, **kw))
    return w


def test_read_writes_once_with_the_outcome_first_and_never_overwrites(tmp_path, ready, capsys):
    out = tmp_path / "reading"
    argv = ["read", "--commit", "abc", "--runs", str(ready.root), "--s3-weight", "plain", "--out", str(out)]
    k30.main(argv)
    body = json.loads((out / "reading.json").read_text())
    assert body["outcome"].startswith("Try K30: NAME dx1") and body["row"] == "K30" and body["tries_counted"] == 3
    assert body["ruling"] == "R53" and body["s3_weight"] == "plain"
    assert body["reading_code"]["script_sha256"] == _sha(T3 / "k30_depth_read.py")
    assert (out / "reading.txt").read_text().splitlines()[0] == body["outcome"]
    shown = capsys.readouterr().out.splitlines()
    assert shown[0] == body["outcome"] and shown[1] == "stop: S1 no, S2 no, S3 no, S4 no"
    assert shown[-1].startswith("wrote ")
    first = {p.name: p.read_bytes() for p in out.iterdir()}
    with pytest.raises(SystemExit) as e:
        k30.main(argv)
    assert e.value.code == 4 and {p.name: p.read_bytes() for p in out.iterdir()} == first
    (out / "reading.json").unlink()  # reading.txt alone is still a reading
    with pytest.raises(SystemExit) as e:
        k30.main(argv)
    assert e.value.code == 4 and not (out / "reading.json").exists()


def test_read_refuses_a_weight_the_row_does_not_name(tmp_path, ready, capsys):
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k30.main(["read", "--commit", "abc", "--runs", str(ready.root), "--s3-weight", "0.5", "--out", str(out)])
    assert e.value.code == 2 and "is not the S3 weight row K30 names (plain)" in capsys.readouterr().err
    assert not out.exists()
    with pytest.raises(SystemExit) as e:  # a weight outside the two
        k30.main(["read", "--commit", "abc", "--runs", str(ready.root), "--s3-weight", "0.7", "--out", str(out)])
    assert e.value.code == 2


def test_an_incomplete_folder_writes_nothing_and_exits_3_until_the_cutoff(tmp_path, ready, monkeypatch, capsys):
    (ready.root / "dm" / "g2_seed_4.json").unlink()
    out = tmp_path / "reading"
    argv = ["read", "--commit", "abc", "--runs", str(ready.root), "--s3-weight", "plain", "--out", str(out)]
    with pytest.raises(SystemExit) as e:
        k30.main(argv)
    assert e.value.code == 3 and not out.exists() and "not a reading" in capsys.readouterr().out
    monkeypatch.setattr(k30, "_now", lambda: datetime(2026, 11, 25, 13, 0, 1, tzinfo=UTC))
    k30.main(argv)
    body = json.loads((out / "reading.json").read_text())
    assert body["void"] and body["outcome"].startswith("Try K30: VOID: INCOMPLETE past its cutoff")
    assert capsys.readouterr().out.splitlines()[0] == body["outcome"]


def test_a_crash_in_read_is_masked_and_writes_nothing(tmp_path, ready, monkeypatch, capsys):
    def boom(runs, **kw):
        raise KeyError("dx1 mean 57.1234 against barred_form 55.4321")

    monkeypatch.setattr(k30, "reading", boom)
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k30.main(["read", "--commit", "abc", "--runs", str(ready.root), "--s3-weight", "plain", "--out", str(out)])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "KeyError" in err and "57.1234" not in err and not out.exists()


def _steps(tmp_path) -> list[list[str]]:
    o = str(tmp_path / "o")
    return [["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--s3-weight", "plain",
             "--out", o],
            ["run", "--commit", "abc", "--world", "dm", "--s3-weight", "plain", "--out", o],
            ["guard", "--commit", "abc", "--runs", o],
            ["read", "--commit", "abc", "--runs", o, "--s3-weight", "plain", "--out", str(tmp_path / "reading")]]


@pytest.mark.parametrize("ledger", ["", ROW, R53, LEDGER.replace("`--s3-weight plain`", "`--s3-weight 0.5`")])
def test_every_step_refuses_without_rows_k30_and_r53_and_loads_nothing(tmp_path, monkeypatch, ledger):
    tries = tmp_path / "TRIES.md"
    tries.write_text(ledger)
    monkeypatch.setattr(k30, "TRIES_MD", tries)
    seen = []
    monkeypatch.setattr(k18, "stamp", lambda: seen.append("stamp") or _prov())
    monkeypatch.setattr(k18, "load_barred_for_plain", lambda *a: seen.append("load"))
    monkeypatch.setattr(k18, "barred_spec", lambda *a: seen.append("spec"))
    for argv in _steps(tmp_path):
        with pytest.raises(SystemExit) as e:
            k30.main(argv)
        assert e.value.code == 2, argv
    assert seen == [] and not (tmp_path / "o").exists() and not (tmp_path / "reading").exists()


@pytest.mark.parametrize("stamp,why", [({"tree_clean": False}, "the tree is not clean"),
                                       ({"commit": "def"}, "HEAD is def, not the run commit abc")])
def test_every_step_refuses_an_unclean_tree_or_another_head(tmp_path, monkeypatch, capsys, stamp, why):
    tries = tmp_path / "TRIES.md"
    tries.write_text(LEDGER)
    monkeypatch.setattr(k30, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: {**_prov(), **stamp})
    for argv in _steps(tmp_path):
        with pytest.raises(SystemExit) as e:
            k30.main(argv)
        assert e.value.code == 2 and why in capsys.readouterr().err, argv
    assert not (tmp_path / "o").exists() and not (tmp_path / "reading").exists()


def test_run_and_guard_refuse_at_and_after_the_cutoff_and_write_nothing(tmp_path, monkeypatch, capsys):
    tries = tmp_path / "TRIES.md"
    tries.write_text(LEDGER)
    monkeypatch.setattr(k30, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: _prov())
    seen = []
    monkeypatch.setattr(k30.k22, "frozen_ok", lambda: seen.append("frozen"))
    for now in (k30.CUTOFF, datetime(2026, 11, 26, 0, 0, tzinfo=UTC)):
        monkeypatch.setattr(k30, "_now", lambda now=now: now)
        for argv in _steps(tmp_path)[:3]:
            with pytest.raises(SystemExit) as e:
                k30.main(argv)
            assert e.value.code == 7 and "past the cutoff (2026-11-25 13:00 UTC)" in capsys.readouterr().err, argv
    assert seen == [] and not (tmp_path / "o").exists()
    monkeypatch.setattr(k30, "_now", lambda: datetime(2026, 11, 25, 12, 59, 59, tzinfo=UTC))  # a second before it
    with pytest.raises(SystemExit) as e:
        k30.main(_steps(tmp_path)[2])
    assert e.value.code != 7 and seen == ["frozen"]  # the clock let it on to the next checks


@pytest.fixture
def checker(tmp_path, monkeypatch, records):
    """cmd_guard over a hand-made guard part: the preconditions passed, K18's rows by hand."""
    w = _write_run(tmp_path / "gc", records, score=False)
    for f in ("guard.json", "guard.json.prov.json"):
        (w.root / f).unlink()
    monkeypatch.setattr(k30, "preconditions", lambda commit: ({"id": "K30", "weight": "plain", "ruling": "R53"},
                                                              _prov()))
    monkeypatch.setattr(k30, "guard_refs", lambda refs=None: w.k18rows)
    return w


def test_the_guard_check_writes_guard_json_once(checker, capsys):
    argv = ["guard", "--commit", "abc", "--runs", str(checker.root)]
    k30.main(argv)
    g = json.loads((checker.root / "guard.json").read_text())
    assert g["passes"] and g["failed"] == [] and g["row"] == "K30" and set(g["files_sha256"]) == set(k30.guard_files())
    assert all(v["distance"]["max_abs_diff"] == 0.0 for v in g["cells"].values())
    assert capsys.readouterr().out == "guard: passes\n"
    with pytest.raises(SystemExit) as e:
        k30.main(argv)
    assert e.value.code == 4


def test_the_guard_check_exits_6_when_the_base_is_off_k18s_rows(tmp_path, monkeypatch, records, capsys):
    w = _write_run(tmp_path / "gf", records, gcells={"dm g2": {k30.BASE: 0.2}}, score=False)
    for f in ("guard.json", "guard.json.prov.json"):
        (w.root / f).unlink()
    monkeypatch.setattr(k30, "preconditions", lambda commit: ({"id": "K30", "weight": "plain", "ruling": "R53"},
                                                              _prov()))
    monkeypatch.setattr(k30, "guard_refs", lambda refs=None: w.k18rows)
    with pytest.raises(SystemExit) as e:
        k30.main(["guard", "--commit", "abc", "--runs", str(w.root)])
    g = json.loads((w.root / "guard.json").read_text())
    assert e.value.code == 6 and g["passes"] is False and g["failed"] == ["barred_form"]
    assert not g["cells"]["dm g2"]["passes"] and "FAILS" in capsys.readouterr().out


@pytest.mark.parametrize("spoil,code,why", [("missing", 3, ""), ("commit", 2, "not this commit's and machine's"),
                                            ("pins", 2, "not at its pinned sha256")])
def test_the_guard_check_refuses_or_waits_and_writes_nothing(checker, monkeypatch, capsys, spoil, code, why):
    f = checker.root / "guard" / "dm" / "g1_seed_6.json"
    if spoil == "missing":
        f.unlink()
    elif spoil == "commit":
        body = json.loads(f.read_text())
        body["provenance"]["commit"] = "def"
        f.write_text(json.dumps(body))
    else:
        def off(refs=None):
            raise k30.Refused("K18's committed dm.json is missing or not at its pinned sha256: refused")

        monkeypatch.setattr(k30, "guard_refs", off)
    with pytest.raises(SystemExit) as e:
        k30.main(["guard", "--commit", "abc", "--runs", str(checker.root)])
    assert e.value.code == code and why in capsys.readouterr().err
    assert not (checker.root / "guard.json").exists()


@pytest.fixture
def fake_run(monkeypatch, tmp_path):
    """cmd_run with the preconditions passed, the barred factor by hand and the harness stubbed: the fake harness logs
    each seed's line in the harness's own format (every arm's score) and leaves the fold's depth vectors where K18's
    plain build leaves them."""
    calls = {"names": [], "seeds": []}
    lines = {}
    st = _prov()
    bs = harness._t3_head()

    def run_task(t, cands, seeds, data_dir, log=print):
        calls["names"] = [c.name for c in cands]
        calls["seeds"] += list(seeds)
        V = _vocab()
        vec = k30.depth_vectors(bs, V, k30._K30["table"])
        k18._FOLD["F"] = {"k30": {"ratio": vec, "V": V, "record": k30.depth_record(V, k30._K30["table"], vec)}}
        for sd in seeds:
            lines[sd] = f"{t.key} seed {sd}: floor 50.00, barred_form {60 - sd / 9:.2f}, dx1 {61 + sd / 7:.2f}"
            log(lines[sd])
        return {"task": {"key": t.key}, "rows": [{"seed": sd} for sd in seeds],
                "provenance": {"commit": "abc", "tree_clean": True}}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=np.ones(3), genes=["a", "b", "c"], spec=SPEC, record={})
        return {"file": "f"}

    monkeypatch.setattr(k30, "preconditions", lambda commit: ({"id": "K30", "weight": "plain", "ruling": "R53"}, st))
    monkeypatch.setattr(k30, "guard_refs", lambda refs=None: {})
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(k18, "stamp", lambda: st)
    out = tmp_path / "run"
    yield SimpleNamespace(calls=calls, lines=lines, out=out, st=st)
    k18._BARRED.clear()


def _run(*extra, part="guard", task="g1", seeds="0-1", weight="plain") -> list[str]:
    return ["run", "--commit", "abc", "--part", part, "--world", "plain", "--task", task, "--seeds", seeds,
            "--s3-weight", weight, *extra]


def test_the_guard_part_runs_the_base_only_and_logs_progress_only(fake_run, capsys):
    argv = _run("--out", str(fake_run.out))
    k30.main(argv)
    assert fake_run.calls["names"] == [harness.FLOOR, "barred_form"] and fake_run.calls["seeds"] == [0, 1]
    printed = capsys.readouterr()
    for x in printed.out.splitlines():
        assert re.fullmatch(r"K30 guard plain g1 seed \d: scored|wrote .+/guard/plain_g1/seed_\d\.json: commit abc, "
                            r"clean True, machine m1", x), x
    for line in fake_run.lines.values():
        for score in re.findall(r"\d+\.\d\d", line):
            assert score not in printed.out and score not in printed.err
    body = json.loads((fake_run.out / "guard" / "plain_g1" / "seed_0.json").read_text())
    assert body["part"] == "guard" and body["arms"] == ["barred_form"] and body["task_source"] == k30.SCRIPT
    assert set(body["depth"]["vectors"]) == {"g1"} and body["depth"]["pairing"] == k30.COLUMN
    assert body["s3_weight"] == "plain" and body["counts"] == {"arms": {}, "pre_step_equals_base": {}}
    assert k18.ARMS == k30._K18_ORIGINAL["ARMS"] and not (fake_run.out / "run_started.json").exists()
    k30.main(argv)  # a rerun skips the seeds written
    assert fake_run.calls["seeds"] == [0, 1] and "present, skipped" in capsys.readouterr().out


def test_a_run_refuses_a_weight_the_row_does_not_name(fake_run, capsys):
    with pytest.raises(SystemExit) as e:
        k30.main(_run("--out", str(fake_run.out), weight="0.5"))
    assert e.value.code == 2 and "is not the S3 weight row K30 names (plain)" in capsys.readouterr().err
    assert fake_run.calls["seeds"] == [] and not fake_run.out.exists()
    with pytest.raises(SystemExit) as e:  # without the flag argparse refuses
        k30.main(["run", "--commit", "abc", "--world", "dm", "--out", str(fake_run.out)])
    assert e.value.code == 2


def test_a_scoring_part_refuses_without_a_passing_guard_and_marks_the_run_steps_start(fake_run, capsys):
    out = fake_run.out
    out.mkdir(parents=True)
    argv = _run("--out", str(out), part="score", task="g2", seeds="0")
    with pytest.raises(SystemExit) as e:
        k30.main(argv)
    assert e.value.code == 2 and "the guard part runs and passes first" in capsys.readouterr().err
    for body, why in (({"passes": False, "provenance": fake_run.st}, "the guard part did not pass"),
                      ({"passes": True, "provenance": {**fake_run.st, "commit": "def"}}, "another commit or machine"),
                      ({"passes": True, "provenance": {**fake_run.st, "machine": "m2"}}, "another commit or machine"),
                      ("[]", "did not pass"), ("not json", "is not JSON")):
        (out / "guard.json").write_text(body if isinstance(body, str) else json.dumps(body))
        with pytest.raises(SystemExit) as e:
            k30.main(argv)
        assert e.value.code == 2 and why in capsys.readouterr().err
    assert fake_run.calls["seeds"] == [] and not (out / "run_started.json").exists()  # no count while refused
    (out / "guard.json").write_text(json.dumps({"passes": True, "provenance": fake_run.st}))
    k30.main(argv)
    assert fake_run.calls["names"] == [harness.FLOOR, *k30.ARMS] and fake_run.calls["seeds"] == [0]
    body = json.loads((out / "plain_g2" / "seed_0.json").read_text())
    assert body["part"] == "score" and body["arms"] == list(k30.ARMS) and set(body["depth"]["vectors"]) == {"g2"}
    mark = json.loads((out / "run_started.json").read_text())
    assert mark["counted"] == ["dx1", "dxh", "dxg"] and mark["commit"] == "abc" and mark["s3_weight"] == "plain"
    assert mark["ruling"] == "R53"


def test_a_table_off_its_pin_refuses_before_any_count(fake_run, synthetic_table, capsys):
    out = fake_run.out
    out.mkdir(parents=True)
    (out / "guard.json").write_text(json.dumps({"passes": True, "provenance": fake_run.st}))
    synthetic_table.write_text(_table_text() + "\n")
    for part in ("guard", "score"):
        with pytest.raises(SystemExit) as e:
            k30.main(_run("--out", str(out), part=part, seeds="0"))
        assert e.value.code == 2 and "not at its pinned sha256" in capsys.readouterr().err
    assert fake_run.calls["seeds"] == [] and not (out / "run_started.json").exists()


def test_a_run_refuses_once_a_reading_is_there_whatever_the_folder(fake_run, monkeypatch, capsys, tmp_path):
    results = tmp_path / "results"
    results.mkdir()
    monkeypatch.setattr(k30, "RESULTS", results)
    argvs = [_run("--out", str(fake_run.out), seeds="0"),
             ["run", "--commit", "abc", "--world", "dm", "--s3-weight", "plain", "--out", str(tmp_path / "another")]]
    for name in ("reading.json", "reading.txt"):
        (results / name).write_text("{}")
        for argv in argvs:
            with pytest.raises(SystemExit) as e:
                k30.main(argv)
            assert e.value.code == 4 and "no arm is scored again" in capsys.readouterr().err, argv
        (results / name).unlink()
    assert fake_run.calls["seeds"] == [] and not fake_run.out.exists() and not (tmp_path / "another").exists()


def test_a_crash_is_masked_writes_no_task_file_and_restores_k18(fake_run, monkeypatch, capsys):
    def boom(*a, **k):
        raise RuntimeError("barred_form 61.234 against dx1 61.2341")

    monkeypatch.setattr(harness, "run_task", boom)
    with pytest.raises(SystemExit) as e:
        k30.main(_run("--out", str(fake_run.out), seeds="0"))
    err = capsys.readouterr().err
    assert e.value.code == 1 and "61.234" not in err and "<num>" in err and "RuntimeError" in err
    assert (k18.ARMS, k18.build_arm, k18.fold_factors) == (k30._K18_ORIGINAL["ARMS"], k30._K18_ORIGINAL["build_arm"],
                                                           k30._K18_ORIGINAL["fold_factors"])
    assert not (fake_run.out / "guard" / "plain_g1" / "seed_0.json").exists()


def test_a_task_flag_on_the_dm_world_or_none_on_the_plain_is_refused():
    for argv in (["run", "--commit", "abc", "--world", "dm", "--task", "g1", "--s3-weight", "plain", "--out", "x"],
                 ["run", "--commit", "abc", "--world", "plain", "--s3-weight", "plain", "--out", "x"]):
        with pytest.raises(SystemExit) as e:
            k30.main(argv)
        assert e.value.code == 2
