#!/usr/bin/env python
"""Try H30 of experiments/harness/TRIES.md (the moonshot specified as "H29: late-anchor chord on the heart panel"; the
label H30 is the ledger's, the label H29 having gone to another round): the heart's expression of kept cells moved
toward a legal late anchor by a chord formed within each outside source, with a target-free platform gate first.

The chord. For each moved lineage (cardiomyocytes; endocardium and endothelium pooled; proepicardium and epicardium,
the lineage map's EpicardiumST; hepatocytes) and each outside source, on the panel genes the source carries:
D = log(mean CP10k at the late anchor + 0.5) - log(mean CP10k at the early anchor + 0.5). The late anchor is JAX
(Qiu 2024, X003) embryos 54-56 and Feng 2022 (X002) E14.5; the early anchor JAX embryo 22 and Feng E9.5 on the one-day
backtest, JAX embryos 13 and 14 alone on the 1.75-day backtest. Genes below 0.5 CP10k at both anchors of a source get
D = 0 in it; the sources' chords are averaged gene by gene where both hold both anchors; |D| is clipped at log 3. A
lineage is moved only with at least 200 kept outside cells at each anchor. Outside cells take their lineage from a
multinomial logistic regression trained on released MERFISH labels (configs/lineage_map.yaml groups) on log1p CP10k
of the panel genes the source carries, with a reject option (top probability below 0.7 or panel library below 30
counts). No author label, cluster or embedding of either source is read: JAX's cell table is read for its embryo_id
and day columns only, Feng's metadata for its stage column only.

Steps (one per call, run_h30.sh):
  gate     STEP G, target-free (it compares no predictor with a fold target and counts no try): JAX's own E8.75 to
           E9.5 change (embryos 13 and 14 to 22) against MERFISH's (E8.75 to E9.5), lineage by lineage, Spearman rho
           over the genes at 1 CP10k or more on both platforms at both stages, against a null of 1000 gene
           permutations within ten bins of MERFISH E8.75's mean log level. PASS, FAIL or VOID; only PASS lets
           anything further run. Writes gate.json and gate.txt; prints one step line.
  tables   STEP K's two chord tables, one per backtest, written with their sha256; refused unless the gate's reading
           of this commit says PASS. Reads no MERFISH E10.5 and no stage after the backtest's reference.
  run      one backtest, a seed range, one JSON: per form X of row H27 (A whole, B half trim, C full trim) the control
           (row H27's form at its forecast size f, then row H25's half), CHORD0 (all factors 1), A (the practice tau),
           A2 (tau 0.4) and P1 to P3 (A2's factors deranged within ten level bins), every arm by
           vec.core.apply_factors(keep_depth=True, seed=[build seed, 29]) on the drawn counts before the size step.
           Refused before anything loads unless rows H30 and N56 are in TRIES.md (row H30 naming the form in force),
           the tree is clean, EXTERNAL_DATA.md's X002 and X003 carry their required wording, the committed gate reading
           says PASS and was made from a clean tree at an ancestor of this commit, the committed table checksums exist
           and the output does not exist. Gate and tables refuse on the same rows, tree and register. Prints no score.
  read     the void conditions, the reproduction gate against row H27's committed scores, the fixed statements on the
           form in force (which the row names) and on the other forms; prints the outcome and what follows only.

The harness's external-window guard (`vec.harness.check_externals`) is bypassed by name: the E9.5 and E8.75 outside
cells lie inside each backtest's analogue window (reference - 1, target + 1]. Every part records it. The read is a
labelled in-sample diagnostic on one truth (E10.5), as T62 and T63 were.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
import time
import zlib
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy import sparse

from vec import harness, prep, provenance, t2_build
from vec.adoption import SEEDS
from vec.boards import GROUPS, load_boards
from vec.core import apply_factors, counts_of, mean_cp10k
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, dense, read_stage
from vec.prep import load_stage, log1p_cp10k

HERE = Path(__file__).resolve().parents[2]  # this checkout
sys.path.insert(0, str(HERE / "experiments" / "external"))


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h28 = _load("h28_heart_trend_strength", HERE / "experiments" / "harness" / "h28_heart_trend_strength.py")
h27, h26, h25, payoff = h28.h27, h28.h26, h28.h25, h28.payoff
import feng2022_geo  # noqa: E402  (the rdata helpers: dgc_matrix, deref, r_class)
import feng_relabel  # noqa: E402  (read_obs: a whitelisted metadata reader; unique_index: symbol -> column)
import filterio as fio  # noqa: E402  (verify_manifest_output)
import windows  # noqa: E402  (parse_stage, allowed_label)

ROW = "H30"
GATE_ROW = "N56"  # the gate's own row (STEP G), committed with row H30 before anything runs
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
REGISTER = HERE / "EXTERNAL_DATA.md"
REGISTER_ROWS = ("X002", "X003")  # each must carry its T2-heart amendment for row H30 before any outside file is read
# the register's required wording, each opening a sentence or a cell of the source's line
AMENDMENT = "T2-heart amendment for row H30"
LICENCE = "Licence confirmed for row H30"
REGISTER_WORDS = {"X002": (AMENDMENT,), "X003": (AMENDMENT, LICENCE)}
AXIS = "h30-heart-late-chord"
PART = "h30.json"
SOURCE = "experiments/harness/h30_heart_late_chord.py (try H30)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h30_heart_late_chord"
COMMITTED_GATE = RESULTS / "gate" / "gate.json"  # the gate's reading, committed with the table checksums
COMMITTED_SHA = RESULTS / "tables_sha256.json"
EXTERNAL = HERE / "data" / "external"
JAX_DIR = "jax_qiu2024"
FENG_DIR = "feng2022_gse193346_geo"

# ------------------------------------------------------------------------------------------------ fixed numbers
SEED = 20261122  # classifier, hold-out split, gate null, split halves; placebo k uses SEED + k
LINEAGES = ("CM", "Endo", "Epi", "Hep")
LINEAGE_OF_GROUP = {"CM": "CM", "Endocardium": "Endo", "Endothelium": "Endo", "EpicardiumST": "Epi",
                    "Hepatocyte": "Hep"}
C_REG, MAX_ITER = 1.0, 2000
HOLDOUT = 0.2
REJECT_PROB = 0.7  # a cell is kept at a top probability of 0.7 or more
MIN_LIBRARY = 30  # and a panel library of 30 counts or more
CHUNK = 50_000  # cells per predict_proba call (memory only; the result does not depend on it)
MIN_CELLS = 200
PSEUDO = 0.5  # log(mean CP10k + 0.5)
GATE_FLOOR = 1.0  # mean CP10k on both platforms at both stages
N_PERM, N_BINS = 1000, 10
BA_MIN, MIN_GENES = 0.80, 100
CM_RHO, OTHER_RHO = 0.35, 0.25
CLIP = float(np.log(3.0))
CHORD_FLOOR = 0.5  # genes below 0.5 CP10k at both anchors of a source get D = 0
TEST_TAU = 0.4
JAX_GATE_EARLY, JAX_GATE_LATE, JAX_GATE_MID = (13, 14), (22,), tuple(range(15, 22))
JAX_LATE = (54, 55, 56)
FENG_LATE = "E14.5"
# the day each JAX embryo must carry in the cell table (the specification's legality check; a range for 15-21)
JAX_DAYS = {13: (8.75, 8.75), 14: (8.75, 8.75), **{n: (9.0, 9.25) for n in range(15, 22)}, 22: (9.5, 9.5),
            54: (14.25, 14.25), 55: (14.333, 14.333), 56: (14.75, 14.75)}
DAY_TOL = 1e-3
TRANSITIONS = {
    "T2h/<=E9.5->E10.5": {"short": "e95", "merfish": "T2/E9.5.h5ad", "jax_early": (22,), "feng_early": "E9.5",
                          "tau": 1.0 / 5.0},
    "T2h/<=E8.75->E10.5": {"short": "e875", "merfish": "T2/E8.75.h5ad", "jax_early": (13, 14), "feng_early": None,
                           "tau": 1.75 / 5.75},
}
TASK_KEYS = tuple(TRANSITIONS)
assert TASK_KEYS == h26.HARNESS_KEYS
# JAX embryos and Feng stages each backtest may never read: inside its forecast span (code guard)
REFUSED_JAX = {"T2h/<=E9.5->E10.5": frozenset(range(23, 54)),
               "T2h/<=E8.75->E10.5": frozenset(range(15, 54))}
REFUSED_FENG = {"T2h/<=E9.5->E10.5": frozenset(), "T2h/<=E8.75->E10.5": frozenset({"E9.5"})}
MERFISH_E105 = "T2/E10.5.h5ad"
GATE_MERFISH = ("T2/E8.75.h5ad", "T2/E9.5.h5ad")

# ------------------------------------------------------------------------------------------------ the arms
FORMS = h27.FORMS  # full_iso_half (A), half_iso_half (B), core_iso_half (C)
LETTER = h26.LETTER
SIZING = "f"
CONTROL, CHORD0, ARM_A, ARM_A2 = "control", "chord0", "A", "A2"
PLACEBOS = ("P1", "P2", "P3")
COUNTED = (ARM_A, ARM_A2, *PLACEBOS)
BUILT = (CONTROL, CHORD0, *COUNTED)
APPLY_TAG = 29  # apply_factors(seed=[build seed, 29])
OWN, CLEAN = h26.OWN, h26.CLEAN
TRUTHS = (OWN, *CLEAN)
WORLDS = (("own", (OWN,)), ("CLEAN", CLEAN))
GUARD_MEAN, GUARD_SEED = h25.GUARD_MEAN, h25.GUARD_SEED  # 0.1 and 0.5 board points
TOL = 1e-9
X_TOL = 1e-5  # CHORD0's largest X difference from the control (row H28 recorded 9.5e-7 for the same pattern)
SHAPE_TERMS = h28.SHAPE_TERMS  # d2_shape, occupancy_dice, scale_log_ratio
DEPTH_TOL = 0.01  # mean per-cell total count change
S1_FOLD, S1_EACH, S2_BAR, MIN_AHEAD = 1.0, 0.5, 0.5, 9
DE_TERMS = ("de_score", "de_direction")
TERMS = ("de_score", "de_direction", "mmd_u", "variogram", "neighborhood_mmd")
DE_GROUPS = ("expression_change", "cell_state")
GROUP_WEIGHT = {name: w for name, w, _ in GROUPS["T2"]}
TAG = "[practice diagnostic, released stages, one truth, in-sample]"
WINDOW_BYPASS = ("vec.harness.check_externals skipped by name: the E9.5 or E8.75 outside cells lie inside the "
                 "backtest's analogue window (reference - 1, target + 1]")


def name(form: str, arm: str) -> str:
    return f"{LETTER[form]}/{arm}"


ARM_NAMES = tuple(name(f, a) for f in FORMS for a in BUILT)
ALL = (harness.FLOOR, *ARM_NAMES)


# ------------------------------------------------------------------------------------------------ refusals
def tries_rows(path: Path = TRIES) -> set[str]:
    return h27.tries_rows(path)


def row_text(row: str = ROW, path: Path = TRIES) -> str | None:
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        if ln.startswith(f"| {row} ") or ln.startswith(f"| {row}|"):
            return ln
    return None


def _opens_sentence(line: str, words: str) -> bool:
    """`words` opens a sentence or a cell of `line` (a phrase quoted inside wording that denies it does not count)."""
    return re.search(rf"(^|[|.;:]\s*){re.escape(words)}(?![\w-])", line) is not None


def register_names_row(path: Path = REGISTER, ids=REGISTER_ROWS, row: str = ROW) -> list[str]:
    """The register rows among `ids` that do not yet carry their required wording for row H30: on X002 and X003
    "T2-heart amendment for row H30", on X003 also "Licence confirmed for row H30", each opening a sentence or a
    cell of the source's line."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    missing = []
    for x in ids:
        ln = next((ln for ln in lines if ln.startswith(f"| {x} ")), None)
        words = tuple(w.replace(ROW, row) for w in REGISTER_WORDS.get(x, (AMENDMENT,)))
        if ln is None or not all(_opens_sentence(ln, w) for w in words):
            missing.append(x)
    return missing


def _clean_or_exit() -> str:
    commit, unclean = provenance.git_state(HERE)
    if unclean is None or unclean:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    return commit


def refuse_rows() -> None:
    """Row H30 and the gate's row N56 both in TRIES.md, and row H30 naming the form in force (TRIES.md cannot change
    between the row commit and the read, so a placeholder there would refuse the read for ever)."""
    have = tries_rows()
    if ROW not in have:
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    if GATE_ROW not in have:
        sys.exit(f"{ROW} is refused: the gate's row {GATE_ROW} is not in experiments/harness/TRIES.md (it is committed "
                 f"with row {ROW} before anything runs)")
    if form_in_force(row_text()) is None:
        sys.exit(f"{ROW} is refused: row {ROW} does not name the form in force ('named in this row before the read: "
                 "<A, B or C>')")


def refuse_common(out_files) -> str:
    """Rows first (H30 naming its form, and the gate's N56), then a clean tree, then the register's amendments, then
    no existing output."""
    refuse_rows()
    commit = _clean_or_exit()
    missing = register_names_row()
    if missing:
        sys.exit(f"{ROW} is refused: EXTERNAL_DATA.md rows {missing} do not carry the required wording for row {ROW} "
                 f"('{AMENDMENT}'; on X003 also '{LICENCE}'), which comes before any outside file is read")
    for p in out_files:
        if Path(p).exists():
            sys.exit(f"{ROW} is refused: {p} exists")
    return commit


def need_rdata() -> None:
    if importlib.util.find_spec("rdata") is None:
        sys.exit(f"{ROW} is refused: the rdata reader is not installed in this environment (uv group external)")


def gate_passed(path: Path, commit: str | None = None) -> bool:
    """The gate's reading at `path` says PASS (and, given a commit, was made at it from a clean tree)."""
    try:
        g = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return False
    if g.get("outcome") != "PASS" or g.get("axis") != AXIS + ":gate":
        return False
    if commit is not None:
        return g.get("provenance", {}).get("commit") == commit and g.get("provenance", {}).get("tree_clean") is True
    return True


def is_ancestor(older: str, newer: str, repo: Path = HERE) -> bool:
    """`older` is a strict ancestor of `newer` (git merge-base --is-ancestor, and not the same commit)."""
    import subprocess

    if not older or not newer or older == newer:
        return False
    r = subprocess.run(["git", "merge-base", "--is-ancestor", older, newer], cwd=repo, capture_output=True)
    return r.returncode == 0


def committed_gate_ok(commit: str) -> bool:
    """The committed gate reading says PASS and was made from a clean tree at a strict ancestor of `commit` (the row
    commit, with the run commit after it)."""
    if not gate_passed(COMMITTED_GATE):
        return False
    try:
        prov = json.loads(Path(COMMITTED_GATE).read_text()).get("provenance", {})
    except (OSError, ValueError):
        return False
    return prov.get("tree_clean") is True and is_ancestor(str(prov.get("commit", "")), str(commit))


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


# ------------------------------------------------------------------------------------------------ outside readers
def _csv_columns(path: Path, columns) -> dict[str, np.ndarray]:
    """The named columns of a CSV and nothing else (each row's other fields are dropped as it is parsed)."""
    with open(path, newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        missing = [c for c in columns if c not in header]
        if missing:
            raise ValueError(f"{Path(path).name}: columns {missing} are not in its header")
        idx = [header.index(c) for c in columns]
        rows = [[r[i] for i in idx] for r in rd]
    arr = np.array(rows, dtype=str).reshape(len(rows), len(idx))
    return {c: arr[:, k] for k, c in enumerate(columns)}


JAX_CELL_COLUMNS = ("embryo_id", "day")  # the only columns of JAX's cell table this script reads


def jax_days(df_cell: Path, embryos) -> dict[int, list[str]]:
    """The day labels of `embryos` in JAX's cell table (its embryo_id and day columns only)."""
    cols = _csv_columns(df_cell, JAX_CELL_COLUMNS)
    want = set(embryos)
    out: dict[int, set[str]] = {n: set() for n in want}
    for e, d in zip(cols["embryo_id"], cols["day"], strict=True):
        m = re.fullmatch(r"embryo_(\d+)", e.strip())
        if m and int(m.group(1)) in want:
            out[int(m.group(1))].add(d.strip())
    return {n: sorted(v) for n, v in out.items()}


def check_jax_days(days: dict[int, list[str]]) -> dict[int, float]:
    """Each embryo carries exactly one day label, the one the specification names (a range for 15-21), outside every
    T2-heart window. Raises otherwise."""
    out = {}
    for n, labels in sorted(days.items()):
        if n not in JAX_DAYS:
            raise ValueError(f"JAX embryo {n} is not one the specification names")
        if len(labels) != 1:
            raise ValueError(f"JAX embryo {n}: day labels {labels}, not exactly one")
        d = windows.parse_stage(labels[0])
        lo, hi = JAX_DAYS[n]
        if not (lo - DAY_TOL <= d <= hi + DAY_TOL):
            raise ValueError(f"JAX embryo {n}: day {labels[0]}, not the specification's E{lo:g}-E{hi:g}")
        if not windows.allowed_label(labels[0], ["T2-heart"]):
            raise ValueError(f"JAX embryo {n}: day {labels[0]} lies inside a T2-heart window")
        out[n] = d
    return out


GENE_ID_COLUMNS = ("gene_ID", "gene_id", "id")
GENE_SYMBOL_COLUMNS = ("gene_short_name", "gene_name", "symbol")


def jax_gene_symbols(df_gene: Path) -> dict[str, str]:
    """Gene ID -> symbol from JAX's gene table (no cells)."""
    with open(df_gene, newline="") as f:
        header = next(csv.reader(f))
    idc = next((c for c in GENE_ID_COLUMNS if c in header), header[0])
    sym = next((c for c in GENE_SYMBOL_COLUMNS if c in header), None)
    if sym is None:
        raise ValueError(f"{Path(df_gene).name}: no symbol column among {GENE_SYMBOL_COLUMNS}")
    cols = _csv_columns(df_gene, (idc, sym))
    return dict(zip(cols[idc], cols[sym], strict=True))


def symbols_of_rows(rows, id_to_symbol: dict[str, str], panel) -> list[str]:
    """The count matrix's row names as symbols: through the gene table when every row name is a gene ID in it, as
    given when the row names already hold the panel's symbols; anything else fails closed."""
    rows = [str(r) for r in rows]
    if rows and all(r in id_to_symbol for r in rows):
        return [id_to_symbol[r] for r in rows]
    if set(panel) & set(rows):
        return rows
    raise ValueError("the count matrix's row names are neither the gene table's IDs nor gene symbols")


def panel_columns(symbols, panel) -> dict[str, int]:
    """Panel gene -> column of a source, for symbols the source holds once (feng_relabel.unique_index)."""
    idx = feng_relabel.unique_index(list(symbols))
    return {g: idx[g] for g in panel if g in idx}


def to_panel(C: sparse.spmatrix, cols: dict[str, int], panel) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Cells x panel genes (integer counts, zero columns for genes the source does not carry) and the carried mask."""
    C = sparse.csr_matrix(C, dtype=np.float64)
    carried = np.array([g in cols for g in panel])
    dest = np.flatnonzero(carried)
    src = np.array([cols[panel[j]] for j in dest], dtype=np.int64)
    S = sparse.csr_matrix((np.ones(len(dest)), (src, dest)), shape=(C.shape[1], len(panel)))
    out = (C @ S).tocsr()
    if out.nnz and not np.allclose(out.data, np.round(out.data)):
        raise ValueError("outside counts are not integers")
    out = sparse.csr_matrix((np.round(out.data).astype(np.int64), out.indices, out.indptr), shape=out.shape)
    out.eliminate_zeros()
    return out, carried


def _parse_rds(path: Path):
    from rdata.parser import parse_file

    return feng2022_geo.deref(parse_file(path).object)


def read_jax_embryo(path: Path, id_to_symbol: dict[str, str], panel) -> tuple[sparse.csr_matrix, np.ndarray]:
    """One JAX embryo's raw counts (a dgCMatrix, genes x cells, in the per-embryo rds) on the panel genes it carries."""
    fio.verify_manifest_output(path)
    obj = _parse_rds(path)
    if feng2022_geo.r_class(obj) != "dgCMatrix":
        raise ValueError(f"{path.name}: expected a dgCMatrix, got {feng2022_geo.r_class(obj)}")
    X, rows, _ = feng2022_geo.dgc_matrix(obj)
    if rows is None:
        raise ValueError(f"{path.name}: the count matrix has no row names")
    cols = panel_columns(symbols_of_rows(rows, id_to_symbol, panel), panel)
    return to_panel(X.T, cols, panel)


def read_jax(root: Path, embryos, panel, refused=frozenset()) -> tuple[dict[int, sparse.csr_matrix], np.ndarray,
                                                                          dict]:
    """JAX embryos' counts on the panel, after the code guard (refused embryos), the embryo-list check against the
    specification's days and the manifest checks."""
    bad = sorted(set(embryos) & set(refused))
    if bad:
        raise PermissionError(f"JAX embryos {bad} lie inside this backtest's forecast span: refused by code")
    root = Path(root)
    df_cell, df_gene = root / "meta" / "df_cell.csv", root / "meta" / "df_gene.csv"
    for p in (df_cell, df_gene):
        fio.verify_manifest_output(p)
    days = check_jax_days(jax_days(df_cell, embryos))
    id_to_symbol = jax_gene_symbols(df_gene)
    out, carried = {}, None
    for n in sorted(embryos):
        C, car = read_jax_embryo(root / "embryo" / f"embryo_{n}_gene_count.rds", id_to_symbol, panel)
        if carried is not None and not np.array_equal(car, carried):
            raise ValueError("JAX embryos carry different panel genes")
        out[n], carried = C, car
    return out, carried, {"embryos": sorted(int(n) for n in embryos), "days": {str(k): v for k, v in days.items()}}


FENG_OBS_COLUMNS = ("stage",)  # the only metadata column of Feng's this script reads


def read_feng(root: Path, stages, panel, refused=frozenset()) -> tuple[dict[str, sparse.csr_matrix], np.ndarray, dict]:
    """Feng's clean CD1 reduction (manifest-checked): the cells of `stages`, raw counts on the panel genes it carries;
    the metadata read for its stage column only (feng_relabel.read_obs, whitelisted)."""
    bad = sorted(set(stages) & set(refused))
    if bad:
        raise PermissionError(f"Feng stages {bad} lie inside this backtest's forecast span: refused by code")
    for s in stages:
        if not windows.allowed_label(s, ["T2-heart"]):
            raise PermissionError(f"Feng stage {s} lies inside a T2-heart window")
    root = Path(root)
    counts, genes_txt, obs = (root / n for n in ("cd1_legal_counts.npz", "cd1_legal_genes.txt", "cd1_legal_obs.csv"))
    for p in (counts, genes_txt, obs):
        fio.verify_manifest_output(p)
    stage = feng_relabel.read_obs(obs, FENG_OBS_COLUMNS)["stage"]
    genes = genes_txt.read_text().split("\n")
    genes = [g for g in genes if g]
    C = sparse.load_npz(counts).tocsr()
    if C.shape != (len(stage), len(genes)):
        raise ValueError("Feng's counts, genes and metadata do not line up")
    cols = panel_columns(genes, panel)
    out = {}
    carried = None
    for s in stages:
        rows = np.flatnonzero(stage == s)
        out[s], carried = to_panel(C[rows], cols, panel)
    return out, carried, {"stages": list(stages)}


# ------------------------------------------------------------------------------------------------ the classifier
def lineage_of(groups) -> np.ndarray:
    return np.array([LINEAGE_OF_GROUP.get(str(g), "") for g in groups], dtype=object).astype(str)


def fit_classifier(C: sparse.spmatrix, labels, seed: int = SEED) -> tuple[object, dict]:
    """sklearn multinomial logistic regression (C 1.0, balanced class weights, max_iter 2000, random_state seed) on
    log1p CP10k of the given genes, trained on 80% of the cells; balanced accuracy on the 20% held out."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.model_selection import train_test_split

    labels = np.asarray(labels).astype(str)
    F = np.asarray(log1p_cp10k(C).todense())
    tr, te = train_test_split(np.arange(len(labels)), test_size=HOLDOUT, random_state=seed, shuffle=True)
    model = LogisticRegression(C=C_REG, class_weight="balanced", max_iter=MAX_ITER, random_state=seed)
    model.fit(F[tr], labels[tr])
    ba = float(balanced_accuracy_score(labels[te], model.predict(F[te])))
    return model, {"balanced_accuracy": ba, "n_train": int(len(tr)), "n_held_out": int(len(te)),
                   "classes": [str(c) for c in model.classes_], "n_iter": int(np.max(model.n_iter_))}


def classify(model, C: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    """Each cell's lineage ('' outside the moved lineages) and its keep mask (top probability 0.7 or more and a panel
    library of 30 counts or more)."""
    C = sparse.csr_matrix(C)
    lib = np.asarray(C.sum(1)).ravel()
    if C.shape[0] == 0:
        return np.array([], dtype=str), np.zeros(0, dtype=bool)
    top, arg = np.empty(C.shape[0]), np.empty(C.shape[0], dtype=np.int64)
    for i in range(0, C.shape[0], CHUNK):  # in chunks: a late whole embryo holds many cells
        P = model.predict_proba(np.asarray(log1p_cp10k(C[i:i + CHUNK]).todense()))
        top[i:i + CHUNK], arg[i:i + CHUNK] = P.max(1), P.argmax(1)
    cls = np.asarray(model.classes_)[arg]
    keep = (top >= REJECT_PROB) & (lib >= MIN_LIBRARY)
    return lineage_of(cls), keep


def log_mean(C: sparse.spmatrix) -> np.ndarray:
    return np.log(mean_cp10k(C) + PSEUDO)


# ------------------------------------------------------------------------------------------------ STEP G
def _bins(basis: np.ndarray, bins: int = N_BINS) -> list[np.ndarray]:
    return np.array_split(np.argsort(np.asarray(basis, float), kind="stable"), bins)


def spearman(a, b) -> float:
    from scipy.stats import spearmanr

    r = spearmanr(np.asarray(a, float), np.asarray(b, float)).statistic
    return float(r) if np.isfinite(r) else float("nan")


def binned_null(x: np.ndarray, y: np.ndarray, basis: np.ndarray, key: str, n: int = N_PERM,
                seed: int = SEED) -> np.ndarray:
    """Spearman rho of x permuted within ten level bins of `basis` against y, n times."""
    rng = np.random.default_rng([seed, zlib.crc32(key.encode())])
    bins = _bins(basis)
    out = np.empty(n)
    for i in range(n):
        p = np.arange(len(x))
        for b in bins:
            p[b] = rng.permutation(b)
        out[i] = spearman(np.asarray(x)[p], y)
    return out


def split_half(early: sparse.spmatrix, late: sparse.spmatrix, genes: np.ndarray, key: str,
               seed: int = SEED) -> float:
    """One platform's split-half ceiling: the delta of one random half of each stage's cells against the other's."""
    rng = np.random.default_rng([seed, zlib.crc32(key.encode())])
    halves = []
    for C in (early, late):
        o = rng.permutation(C.shape[0])
        halves.append((C[o[: len(o) // 2]], C[o[len(o) // 2:]]))
    d1 = log_mean(halves[1][0])[genes] - log_mean(halves[0][0])[genes]
    d2 = log_mean(halves[1][1])[genes] - log_mean(halves[0][1])[genes]
    return spearman(d1, d2)


def gate_lineage(lin: str, m_early, m_late, j_early, j_late, j_mid=None) -> dict:
    """One lineage of STEP G: MERFISH and JAX cells (counts on the genes JAX carries) at E8.75 and E9.5 (JAX: embryos
    13 and 14, embryo 22; j_mid embryos 15-21, reported only)."""
    sets = {"merfish_E8.75": m_early, "merfish_E9.5": m_late, "jax_E8.75": j_early, "jax_E9.5": j_late}
    n = {k: int(v.shape[0]) for k, v in sets.items()}
    rec = {"lineage": lin, "cells": n, "evaluable": all(v >= MIN_CELLS for v in n.values())}
    if not rec["evaluable"]:
        return {**rec, "passes": False}
    means = {k: mean_cp10k(v) for k, v in sets.items()}
    floor = np.all([means[k] >= GATE_FLOOR for k in sets], axis=0)
    genes = np.flatnonzero(floor)
    rec["genes"] = int(len(genes))
    if len(genes) < MIN_GENES:
        return {**rec, "passes": False, "too_few_genes": True}
    dm = np.log(means["merfish_E9.5"][genes] + PSEUDO) - np.log(means["merfish_E8.75"][genes] + PSEUDO)
    dj = np.log(means["jax_E9.5"][genes] + PSEUDO) - np.log(means["jax_E8.75"][genes] + PSEUDO)
    basis = np.log(means["merfish_E8.75"][genes] + PSEUDO)
    rho = spearman(dj, dm)
    null = binned_null(dj, dm, basis, lin)
    p95 = float(np.percentile(null, 95))
    thr = CM_RHO if lin == "CM" else OTHER_RHO
    rec.update(rho=rho, null_p95=p95, null_mean=float(null.mean()), threshold=thr,
               passes=bool(rho >= thr and rho > p95), too_few_genes=False,
               ceiling={"merfish": split_half(m_early, m_late, genes, f"{lin}|merfish"),
                        "jax": split_half(j_early, j_late, genes, f"{lin}|jax")})
    if j_mid is not None and j_mid.shape[0] >= MIN_CELLS:
        dmid = np.log(mean_cp10k(j_mid)[genes] + PSEUDO) - np.log(means["jax_E8.75"][genes] + PSEUDO)
        rec["intermediate"] = {"rho": spearman(dmid, dm), "cells": int(j_mid.shape[0])}
    else:
        rec["intermediate"] = {"rho": None, "cells": int(0 if j_mid is None else j_mid.shape[0])}
    return rec


def gate_outcome(lineages: dict[str, dict], ba: float) -> tuple[str, list[str]]:
    """VOID (balanced accuracy below 0.80, or an evaluable lineage with fewer than 100 genes over the floor), PASS
    (CM passes and at least one other lineage passes) or FAIL."""
    why = []
    if not ba >= BA_MIN:
        why.append(f"the classifier's balanced accuracy on held-out MERFISH is below {BA_MIN}")
    few = [k for k, r in lineages.items() if r.get("too_few_genes")]
    if few:
        why.append(f"fewer than {MIN_GENES} genes pass the floor in {few}")
    if why:
        return "VOID", why
    cm = lineages.get("CM", {}).get("passes", False)
    other = any(r.get("passes") for k, r in lineages.items() if k != "CM")
    return ("PASS" if cm and other else "FAIL"), []


def by_lineage(C: sparse.spmatrix, lin: np.ndarray, keep: np.ndarray | None = None) -> dict[str, sparse.csr_matrix]:
    C = sparse.csr_matrix(C)
    keep = np.ones(C.shape[0], bool) if keep is None else keep
    return {k: C[np.flatnonzero(keep & (lin == k))] for k in LINEAGES}


def merfish_stage(data_dir: Path, f: str, genes) -> tuple[sparse.csr_matrix, np.ndarray]:
    """A released MERFISH stage's counts (duplicates dropped) on `genes` and its lineage-map groups. E10.5 refused."""
    if f == MERFISH_E105:
        raise PermissionError("MERFISH E10.5 is the backtests' target: no table or gate reads it")
    a = load_stage(Path(data_dir) / f, genes=list(genes), drop_duplicates=True)
    return counts_of(a, list(genes)), np.asarray(a.obs["group"]).astype(str)


GATE_FOLLOWS = {
    "PASS": "gate PASS: the platform premise holds by the fixed line; STEP K may run (its tables first).",
    "FAIL": "gate FAIL: the platform premise is dead by the fixed line; H30 ends, no counted read runs, nothing is "
            "asked.",
    "VOID": "gate VOID: {why}; H30 ends, no counted read runs, nothing is asked.",
}


def run_gate(out: Path, data_dir: Path, external: Path) -> dict:
    out = Path(out)
    commit = refuse_common([out / "gate.json", out / "gate.txt"])
    need_rdata()
    boards = load_boards()
    panel = list(boards[payoff.heart_task(TASK_KEYS[0]).panel_board].genes)
    jax, carried, jrec = read_jax(Path(external) / JAX_DIR, (*JAX_GATE_EARLY, *JAX_GATE_MID, *JAX_GATE_LATE), panel)
    genes = [g for g, c in zip(panel, carried, strict=True) if c]
    gi = np.flatnonzero(carried)
    m = {f: merfish_stage(data_dir, f, panel) for f in GATE_MERFISH}
    Ctrain = sparse.vstack([m[f][0][:, gi] for f in GATE_MERFISH]).tocsr()
    model, crec = fit_classifier(Ctrain, np.concatenate([m[f][1] for f in GATE_MERFISH]))

    def jcells(embryos):
        C = sparse.vstack([jax[n][:, gi] for n in embryos]).tocsr()
        lin, keep = classify(model, C)
        return by_lineage(C, lin, keep)

    je, jl, jm = jcells(JAX_GATE_EARLY), jcells(JAX_GATE_LATE), jcells(JAX_GATE_MID)
    me, ml = (by_lineage(m[f][0][:, gi], lineage_of(m[f][1])) for f in GATE_MERFISH)
    lineages = {k: gate_lineage(k, me[k], ml[k], je[k], jl[k], jm[k]) for k in LINEAGES}
    outcome, why = gate_outcome(lineages, crec["balanced_accuracy"])
    res = {"axis": AXIS + ":gate", "row": ROW, "outcome": outcome, "why": why, "classifier": crec,
           "genes_matched": len(genes), "genes": genes, "jax": jrec, "merfish": list(GATE_MERFISH),
           "lineages": lineages, "fixed": {"seed": SEED, "reject_prob": REJECT_PROB, "min_library": MIN_LIBRARY,
                                            "min_cells": MIN_CELLS, "floor_cp10k": GATE_FLOOR, "perm": N_PERM,
                                            "bins": N_BINS, "ba_min": BA_MIN, "min_genes": MIN_GENES,
                                            "cm_rho": CM_RHO, "other_rho": OTHER_RHO},
           "note": "target-free: compares no predictor with a fold target; counts no try",
           "source": SOURCE, "provenance": {**payoff.provenance(), "commit_at_start": commit}}
    res["follows"] = GATE_FOLLOWS[outcome].format(why="; ".join(why))
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "gate.txt").write_text(gate_text(res))
    print(f"gate {outcome}", flush=True)
    return res


def gate_text(res: dict) -> str:
    lines = [f"Try H30, STEP G (commit {res['provenance']['commit']}): target-free, counts no try. "
             f"OUTCOME: {res['outcome']}.", res["follows"],
             f"Classifier balanced accuracy on held-out MERFISH: {res['classifier']['balanced_accuracy']:.3f} "
             f"(void below {BA_MIN}).", f"Panel genes JAX carries: {res['genes_matched']} (the specification: 475)."]
    for k, r in res["lineages"].items():
        if not r["evaluable"]:
            lines.append(f"{k}: not evaluable (cells {r['cells']}; at least {MIN_CELLS} each); does not pass")
        elif r.get("too_few_genes"):
            lines.append(f"{k}: {r['genes']} genes over the floor (fewer than {MIN_GENES}): void")
        else:
            mid = r["intermediate"]["rho"]
            lines.append(f"{k}: rho {r['rho']:+.3f} (line {r['threshold']}, null 95th percentile {r['null_p95']:+.3f})"
                         f" {'PASSES' if r['passes'] else 'does not pass'}; genes {r['genes']}; cells {r['cells']}; "
                         f"split-half ceiling MERFISH {r['ceiling']['merfish']:+.3f}, JAX {r['ceiling']['jax']:+.3f}; "
                         f"embryos 15-21 as the late point (reported only) "
                         f"{'n/a' if mid is None else f'{mid:+.3f}'}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ STEP K tables
def source_chord(early: sparse.spmatrix, late: sparse.spmatrix, carried: np.ndarray) -> np.ndarray:
    """One source's chord on the panel (NaN where it carries no gene): log(mean CP10k late + 0.5) - log(mean CP10k
    early + 0.5); 0 where both anchors are below 0.5 CP10k."""
    me, ml = mean_cp10k(early), mean_cp10k(late)
    d = np.log(ml + PSEUDO) - np.log(me + PSEUDO)
    d[(me < CHORD_FLOOR) & (ml < CHORD_FLOOR)] = 0.0
    d[~np.asarray(carried, bool)] = np.nan
    return d


def combine(chords: list[np.ndarray], n_genes: int) -> np.ndarray:
    """The sources' chords averaged gene by gene where present, 0 where none is, clipped at |log 3|."""
    if not chords:
        return np.zeros(n_genes)
    S = np.vstack(chords)
    have = np.isfinite(S)
    d = np.divide(np.where(have, S, 0.0).sum(0), have.sum(0), out=np.zeros(n_genes), where=have.sum(0) > 0)
    return np.clip(d, -CLIP, CLIP)


def chord_table(anchor_cells: dict[str, dict[str, dict[str, sparse.csr_matrix]]],
                carried: dict[str, np.ndarray], n_genes: int) -> tuple[dict[str, np.ndarray], dict]:
    """anchor_cells[source][anchor]['CM'...]: kept outside cells of each lineage at 'early' and 'late'. A source
    contributes to a lineage where it holds both anchors with at least one kept cell at each; the lineage is moved
    only with at least 200 kept cells at each anchor over the sources."""
    D, rec = {}, {}
    for lin in LINEAGES:
        n = {s: {a: int(anchor_cells[s][a][lin].shape[0]) for a in anchor_cells[s]} for s in anchor_cells}
        total = {a: sum(n[s].get(a, 0) for s in n) for a in ("early", "late")}
        contrib = sorted(s for s in anchor_cells if set(anchor_cells[s]) == {"early", "late"}
                         and n[s]["early"] > 0 and n[s]["late"] > 0)
        moved = bool(contrib) and all(v >= MIN_CELLS for v in total.values())
        rec[lin] = {"cells": n, "total": total, "sources": contrib, "moved": moved}
        if moved:
            D[lin] = combine([source_chord(anchor_cells[s]["early"][lin], anchor_cells[s]["late"][lin], carried[s])
                              for s in contrib], n_genes)
    return D, rec


def build_table(key: str, data_dir: Path, external: Path, panel, jax_cache: dict | None = None) -> dict:
    """One backtest's table: classifier on the reference stage's MERFISH labels only, per source on the panel genes the
    source carries; chord per lineage. Reads no MERFISH E10.5 and no stage after the reference."""
    tr = TRANSITIONS[key]
    t = payoff.heart_task(key)
    if tr["merfish"] != t.reference or tr["merfish"] == MERFISH_E105:
        raise PermissionError(f"{key}: the classifier reads the reference stage only")
    Cm, groups = merfish_stage(data_dir, tr["merfish"], panel)
    jax_emb = (*tr["jax_early"], *JAX_LATE)
    jcache = jax_cache if jax_cache is not None else {}
    bad = sorted(set(jax_emb) & REFUSED_JAX[key])
    if bad:
        raise PermissionError(f"JAX embryos {bad} lie inside {key}'s forecast span: refused by code")
    need = [n for n in jax_emb if n not in jcache.get("counts", {})]
    if need:
        got, car, jrec = read_jax(Path(external) / JAX_DIR, need, panel, REFUSED_JAX[key])
        if "carried" in jcache and not np.array_equal(jcache["carried"], car):
            raise ValueError("JAX embryos carry different panel genes")
        jcache.setdefault("counts", {}).update(got)
        jcache.setdefault("days", {}).update({int(k): v for k, v in jrec["days"].items()})
        jcache["carried"] = car
    days = {n: jcache["days"][n] for n in jax_emb}
    feng_stages = tuple(s for s in (tr["feng_early"], FENG_LATE) if s)
    fc, fcar, _ = read_feng(Path(external) / FENG_DIR, feng_stages, panel, REFUSED_FENG[key])
    carried = {"jax": jcache["carried"], "feng": fcar}
    anchors = {"jax": {"early": tr["jax_early"], "late": JAX_LATE},
               "feng": {**({"early": tr["feng_early"]} if tr["feng_early"] else {}), "late": FENG_LATE}}
    cells, classifiers = {}, {}
    for src in ("jax", "feng"):
        gi = np.flatnonzero(carried[src])
        model, crec = fit_classifier(Cm[:, gi], groups)
        classifiers[src] = {**crec, "genes": int(len(gi)), "trained_on": [tr["merfish"]]}
        cells[src] = {}
        for a, what in anchors[src].items():
            C = (sparse.vstack([jcache["counts"][n] for n in what]).tocsr() if src == "jax" else fc[what])
            lin, keep = classify(model, C[:, gi])
            cells[src][a] = by_lineage(C, lin, keep)
    D, rec = chord_table(cells, carried, len(panel))
    return {"axis": AXIS + ":table", "row": ROW, "task": key, "tau_practice": tr["tau"], "tau_test": TEST_TAU,
            "genes": list(panel), "lineages": {k: [float(x) for x in v] for k, v in D.items()},
            "record": rec, "classifiers": classifiers,
            "anchors": {"jax": {a: [int(n) for n in v] for a, v in anchors["jax"].items()},
                        "feng": anchors["feng"]},
            "jax_days": {str(k): v for k, v in days.items()},
            "carried": {s: int(np.sum(v)) for s, v in carried.items()},
            "stages_read": {"merfish": [tr["merfish"]], "jax": sorted(int(n) for n in jax_emb),
                            "feng": list(feng_stages)},
            "fixed": {"clip": CLIP, "floor": CHORD_FLOOR, "pseudo": PSEUDO, "min_cells": MIN_CELLS, "seed": SEED},
            "source": SOURCE, "provenance": payoff.provenance()}


def table_bytes(table: dict) -> bytes:
    return (json.dumps(table, indent=1, sort_keys=True, default=float) + "\n").encode()


def run_tables(out: Path, gate_dir: Path, data_dir: Path, external: Path) -> dict:
    out = Path(out)
    files = [out / f"{TRANSITIONS[k]['short']}.json" for k in TASK_KEYS] + [out / "tables_sha256.json"]
    commit = refuse_common(files)
    if not gate_passed(Path(gate_dir) / "gate.json", commit):
        sys.exit(f"{ROW} is refused: the gate's reading of this commit does not say PASS")
    need_rdata()
    boards = load_boards()
    panel = list(boards[payoff.heart_task(TASK_KEYS[0]).panel_board].genes)
    out.mkdir(parents=True, exist_ok=True)
    cache: dict = {}
    sums = {}
    for key in TASK_KEYS:
        b = table_bytes(build_table(key, data_dir, external, panel, cache))
        f = out / f"{TRANSITIONS[key]['short']}.json"
        f.write_bytes(b)
        sums[f.name] = sha256_bytes(b)
        print(f"tables {TRANSITIONS[key]['short']} written", flush=True)
    (out / "tables_sha256.json").write_text(json.dumps(sums, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out / 'tables_sha256.json'}", flush=True)
    return sums


def load_table(path: Path, committed: Path | None = None) -> tuple[dict, str]:
    """A table, refused unless its sha256 is the committed one."""
    path = Path(path)
    committed = COMMITTED_SHA if committed is None else committed
    b = path.read_bytes()
    got = sha256_bytes(b)
    try:
        want = json.loads(Path(committed).read_text()).get(path.name)
    except (OSError, ValueError):
        want = None
    if want is None or want != got:
        sys.exit(f"{ROW} is refused: {path.name}'s sha256 is not the committed one")
    return json.loads(b), got


# ------------------------------------------------------------------------------------------------ the arms
def factor_rows(table: dict, genes, tau: float) -> dict[str, np.ndarray]:
    """exp(tau D) per moved lineage, in `genes` order (1 for a gene the table does not hold)."""
    pos = {g: i for i, g in enumerate(table["genes"])}
    idx = np.array([pos.get(g, -1) for g in genes])
    out = {}
    for lin, d in table["lineages"].items():
        d = np.asarray(d, float)
        v = np.where(idx >= 0, d[np.maximum(idx, 0)], 0.0)
        out[lin] = np.exp(tau * v)
    return out


def placebo_rows(rows: dict[str, np.ndarray], basis: dict[str, np.ndarray], k: int) -> tuple[dict, dict]:
    """A2's factors with each lineage's genes deranged within ten bins of its reference-stage MERFISH mean log level
    (row H28's derangement at seed SEED + k)."""
    out, fixed = {}, {}
    for lin, f in rows.items():
        p = h28.derangement(basis[lin], lin, bins=N_BINS, seed=SEED + k)
        out[lin] = np.asarray(f)[p]
        fixed[lin] = int((p == np.arange(len(p))).sum())
    return out, fixed


def reference_basis(ref: prep.ad.AnnData, genes, lineages) -> dict[str, np.ndarray]:
    """log(mean CP10k + 0.5) of each lineage's cells (released labels) in the whole reference stage."""
    C = counts_of(ref, list(genes))
    lin = lineage_of(np.asarray(ref.obs["group"]).astype(str))
    return {k: log_mean(C[np.flatnonzero(lin == k)]) if (lin == k).any() else np.zeros(len(genes)) for k in lineages}


def chord_expression(ctx: harness.BuildContext, rows_by_lineage: dict[str, np.ndarray]):
    """The context's draw with count-space factors per moved lineage (released labels of the drawn cells) through
    apply_factors(keep_depth=True, seed=[build seed, 29]); X re-rendered as log1p CP10k of the new counts, then
    harness.finish. Returns the stage, the new counts, the draw's rows and a record."""
    pool, rows, _, coords = harness.draw_copy(ctx)
    genes = list(ctx.board.genes)
    drawn = counts_of(pool[rows], genes)
    lin = lineage_of(np.asarray(pool.obs["group"])[rows])
    moved = sorted(rows_by_lineage)
    F = np.vstack([np.ones(len(genes)), *(rows_by_lineage[k] for k in moved)])
    ti = np.zeros(len(rows), dtype=np.int64)
    for i, k in enumerate(moved):
        ti[lin == k] = i + 1
    new, scale, dead = apply_factors(drawn, F, ti, keep_depth=True, seed=[ctx.seed, APPLY_TAG])
    st = harness.finish(ctx, pool, dense(log1p_cp10k(new)), coords)
    t0, t1 = np.asarray(drawn.sum(1)).ravel().astype(float), np.asarray(new.sum(1)).ravel().astype(float)
    rec = {"moved": moved, "moved_cells": {k: int((ti == i + 1).sum()) for i, k in enumerate(moved)},
           "mean_total_change": float(np.mean(np.divide(t1 - t0, t0, out=np.zeros_like(t0), where=t0 > 0))),
           "restored_cells": int(dead.sum()), "unmoved_counts_equal": bool((new[ti == 0] != drawn[ti == 0]).nnz == 0)}
    return st, new, rows, rec


def build_arms(ctx: harness.BuildContext, form: str, size: dict, factors: dict[str, dict[str, np.ndarray]],
               fixed: dict) -> tuple[dict, dict]:
    """One form's arms on its pool's context: the control (row H27's form at f), CHORD0, A, A2, P1 to P3, each sized
    by row H27's f step, halved and checked; and the identity checks."""
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    pre = {CONTROL: (harness._copy(ctx), None, None, {})}
    for arm in (CHORD0, *COUNTED):
        pre[arm] = chord_expression(ctx, factors[arm])
    pool, rows, _, _ = harness.draw_copy(ctx)
    drawn = counts_of(pool[rows], list(ctx.board.genes))
    ctrl = pre[CONTROL][0]
    checks = {"chord0_counts_equal_control": bool((pre[CHORD0][1] != drawn).nnz == 0),
              "same_draw": bool(all(np.array_equal(v[2], rows) for k, v in pre.items() if k != CONTROL)),
              "pre_size_coords_equal": bool(all(np.array_equal(v[0].coords, ctrl.coords) for v in pre.values())),
              "chord0_x_max_abs_diff_control": float(np.abs(np.asarray(pre[CHORD0][0].X, np.float64)
                                                            - np.asarray(ctrl.X, np.float64)).max()),
              "placebo_fixed_points": fixed}
    out, recs = {}, {}
    for arm, (st, _, _, rec) in pre.items():
        X, C, srec = h27.size_arm(st.X, st.coords, form, SIZING, {}, size, {}, seed)
        sub, rep = t2_build.build_submission(ctx.board, X, C, genes=list(ctx.board.genes))
        if not rep.ok:
            raise ValueError(f"{ctx.task.key} {form} {arm}: the built file fails its check: {rep.errors}")
        out[arm] = Stage.from_anndata(sub, ctx.board.genes, need_coords=True)
        recs[arm] = {"size": srec, "chord": rec}
    checks["post_size_coords_equal"] = bool(all(np.array_equal(v.coords, out[CONTROL].coords) for v in out.values()))
    checks["n_cells"] = {k: int(v.n) for k, v in out.items()}
    return out, {"records": recs, "checks": checks}


def arm_factors(table: dict, genes, basis: dict[str, np.ndarray]) -> tuple[dict, dict]:
    """Every arm's per-lineage factor rows: CHORD0 ones, A at the practice tau, A2 at 0.4, P1-P3 A2 deranged."""
    a2 = factor_rows(table, genes, TEST_TAU)
    out = {CHORD0: {k: np.ones(len(genes)) for k in a2}, ARM_A: factor_rows(table, genes, table["tau_practice"]),
           ARM_A2: a2}
    fixed = {}
    for k, p in enumerate(PLACEBOS, start=1):
        out[p], fixed[p] = placebo_rows(a2, basis, k)
    return out, fixed


# ------------------------------------------------------------------------------------------------ the run
def refuse_run(path: Path) -> str:
    commit = refuse_common([path])
    if not gate_passed(COMMITTED_GATE):
        sys.exit(f"{ROW} is refused: the committed gate reading does not say PASS (FAIL or VOID ends H30)")
    if not committed_gate_ok(commit):
        sys.exit(f"{ROW} is refused: the committed gate reading was not made from a clean tree at an ancestor of this "
                 "commit")
    if not COMMITTED_SHA.exists():
        sys.exit(f"{ROW} is refused: the tables' sha256 is not committed")
    return commit


def run_transition(key: str, seeds, out: Path, data_dir: Path, tables: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse_run(path)
    if key not in TASK_KEYS:
        sys.exit(f"{key}: not a backtest of row {ROW}")
    table, table_sha = load_table(Path(tables) / f"{TRANSITIONS[key]['short']}.json")
    if table.get("task") != key:
        sys.exit(f"{ROW} is refused: the table is not {key}'s")
    t = payoff.heart_task(key)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    geoms, _ = h27.input_geometries(t, data_dir)  # the forecast first, from the input stages only (row H27's)
    size = h27.analogue_size_result(geoms, float(harness.STAGE_DAYS[t.target]))
    ftargets = h27.forecast_targets(size)
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    tm = h26.truth_masks(TC, tcore)
    truths = {q: Stage.from_anndata(ta[tm[q]], genes, need_coords=True) for q in TRUTHS}
    t_info = {q: h26.truth_info(TC, tcore, tm[q]) for q in TRUTHS}
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    whole = inputs["ref"]
    pcore = np.asarray(whole.obs["core"], dtype=bool)
    PC = np.asarray(whole.obsm[COORD_KEY], dtype=float)[:, :3]
    pm = {f: m for f, m in h26.pool_masks(PC, pcore).items() if f in FORMS}
    pools = {f: (whole if h26.TRIM[f] is None else whole[m].copy()) for f, m in pm.items()}
    p_info = {f: h26.pool_info(PC, pcore, m) for f, m in pm.items()}
    v7 = h27.v7_crosscheck(PC, size)
    basis = reference_basis(whole, list(genes), sorted(table["lineages"]))
    factors, fixed = arm_factors(table, list(genes), basis)
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        preds, recs, checks, n_cells = {}, {}, {}, {}
        for form in FORMS:
            ctx = harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pools[form]},
                                       seed=harness.BUILD_SEED_OFFSET + s)
            with harness.guarded_reads(harness.REPO_ROOT, data_dir):
                built, extra = build_arms(ctx, form, size, factors, fixed)
            for arm, st in built.items():
                preds[name(form, arm)] = st
                recs[name(form, arm)] = extra["records"][arm]
                n_cells[name(form, arm)] = extra["checks"]["n_cells"][arm]
            checks[LETTER[form]] = extra["checks"]
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": n_cells, "checks": checks,
               "records": recs, "scores": {}, "skills": {}, "groups": {}}
        for q, truth in truths.items():
            board = replace(base_board, truth_cells=truth.n)
            fl = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, s), truth, ref_stage, seed=s,
                             anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
            row["scores"][q] = {harness.FLOOR: fl.score}
            row["skills"][q] = {harness.FLOOR: fl.skills}
            row["groups"][q] = {harness.FLOOR: dict(fl.group_skill)}
            for nm, pred in preds.items():
                r = score_board(board, pred, truth, ref_stage, seed=s, anchors=fl.anchors, build=harness.BUILD,
                                protocol=harness.PROTOCOL)
                row["scores"][q][nm], row["skills"][q][nm] = r.score, r.skills
                row["groups"][q][nm] = dict(r.group_skill)
        row["seconds"] = round(time.time() - t0, 1)
        rows_out.append(row)
        print(f"{key} seed {s} ({row['seconds']} s) done", flush=True)
    res = {"label": "try H30: the late-anchor chord on the heart panel, three forms; a labelled in-sample diagnostic "
                    "on one truth, not an adoption result", "axis": AXIS, "task": t.key, "arms": list(ARM_NAMES),
           "truth_values": list(TRUTHS), "forecast_targets": ftargets, "size_result": size, "v7": v7,
           "truth_info": t_info, "pred_info": p_info, "truth_seed": payoff.TRUTH_SEED,
           "half_seed_base": h25.HALF_SEED_BASE, "table_sha256": table_sha,
           "table": {k: table[k] for k in ("task", "tau_practice", "tau_test", "record", "classifiers", "anchors",
                                           "jax_days", "stages_read", "carried")},
           "external_window": list(harness.external_window(t)), "window_guard": WINDOW_BYPASS,
           "seeds": list(seeds), "rows": rows_out, "source": SOURCE, "provenance": payoff.provenance()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


# ------------------------------------------------------------------------------------------------ loading
def load_parts(root: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for p in sorted(Path(root).rglob(PART)):
        r = json.loads(p.read_text())
        if r.get("axis") == AXIS:
            out.setdefault(r["task"], []).append(r)
    return out


def load_h27(root: Path = h28.H27_PARTS) -> dict[str, list[dict]]:
    return h27.load_parts(root)


# ------------------------------------------------------------------------------------------------ the void
def _expected_table(key: str) -> dict:
    tr = TRANSITIONS[key]
    return {"merfish": [tr["merfish"]], "jax": sorted((*tr["jax_early"], *JAX_LATE)),
            "feng": [s for s in (tr["feng_early"], FENG_LATE) if s]}


def table_checks(key: str, part: dict, sha: dict) -> list[str]:
    """V6: the table's sha256 is the committed one; it read the reference stage only (no MERFISH E10.5, no stage
    between the reference and the target), exactly the specification's JAX embryos and Feng stages, at the
    specification's taus."""
    why = []
    tr = TRANSITIONS[key]
    if sha.get(f"{tr['short']}.json") != part.get("table_sha256"):
        why.append(f"V6: {key}: the table's sha256 is not the committed one")
    tb = part["table"]
    got = tb["stages_read"]
    want = _expected_table(key)
    if got.get("merfish") != want["merfish"]:
        why.append(f"V6: {key}: a fit read MERFISH {got.get('merfish')}, not the reference stage alone")
    ref_day, tgt = harness.STAGE_DAYS[tr["merfish"]], harness.STAGE_DAYS[MERFISH_E105]
    for f in [*got.get("merfish", []), *(c for s in tb["classifiers"].values() for c in s["trained_on"])]:
        d = harness.STAGE_DAYS.get(f)
        if f == MERFISH_E105 or d is None or ref_day < d <= tgt:
            why.append(f"V6: {key}: a fit read {f}")
    if sorted(got.get("jax", [])) != want["jax"]:
        why.append(f"V6: {key}: the JAX embryo list {got.get('jax')} is not {want['jax']}")
    if list(got.get("feng", [])) != want["feng"]:
        why.append(f"V6: {key}: the Feng stages {got.get('feng')} are not {want['feng']}")
    if not (abs(tb["tau_practice"] - tr["tau"]) <= 1e-15 and tb["tau_test"] == TEST_TAU):
        why.append(f"V6: {key}: the taus are not the specification's")
    return why


def identities(key: str, rows: list[dict]) -> list[str]:
    """V5: CHORD0's counts are the control's and its X within 1e-5 of the control's; every arm has the control's draw, coordinates (before and after the size
    step and half) and cell number; no shape term moves by more than 1e-9; no placebo gene keeps its factor; the mean
    per-cell total count change is within 1%; unmoved cells keep their counts."""
    why = []
    for r in rows:
        tag = f"V5: {key} seed {r['seed']}"
        for form in FORMS:
            c = r["checks"][LETTER[form]]
            if c.get("chord0_counts_equal_control") is not True:
                why.append(f"{tag} {LETTER[form]}: CHORD0's counts are not the control's")
            xd = c.get("chord0_x_max_abs_diff_control")
            if isinstance(xd, bool) or not isinstance(xd, (int, float)) or not abs(xd) <= X_TOL:
                why.append(f"{tag} {LETTER[form]}: CHORD0's X differs from the control's by more than {X_TOL:g}")
            for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
                if c.get(k) is not True:
                    why.append(f"{tag} {LETTER[form]}: {k.replace('_', ' ')} fails")
            if any(v for p in c.get("placebo_fixed_points", {}).values() for v in p.values()):
                why.append(f"{tag} {LETTER[form]}: a placebo leaves a gene on its own factor")
            names = [name(form, a) for a in BUILT]
            if len({r["n_cells"][n] for n in names}) != 1:
                why.append(f"{tag} {LETTER[form]}: the arms' cell numbers differ")
            for a in (CHORD0, *COUNTED):
                rec = r["records"][name(form, a)]["chord"]
                if not abs(rec["mean_total_change"]) <= DEPTH_TOL:
                    why.append(f"{tag} {name(form, a)}: mean per-cell total count change above 1%")
                if rec.get("unmoved_counts_equal") is not True:
                    why.append(f"{tag} {name(form, a)}: an unmoved cell's counts changed")
            for q in TRUTHS:
                for a in BUILT:
                    for m in SHAPE_TERMS:
                        if abs(r["skills"][q][name(form, a)][m] - r["skills"][q][name(form, CONTROL)][m]) > TOL:
                            why.append(f"{tag} {q}: {name(form, a)}'s {m} is not the control's")
    return why


def void_reasons(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], sha: dict,
                 commit: str | None = None) -> tuple[list[str], dict]:
    why: list[str] = []
    dist: dict = {"V4": {}}
    every = [p for ps in parts.values() for p in ps]
    complete = True
    for key in TASK_KEYS:  # V1
        ps = parts.get(key)
        if not ps:
            why.append(f"V1: missing transition {key}")
            complete = False
            continue
        rows, twice = h26.merge_rows(ps)
        if twice:
            why.append(f"V1: {key}: seeds {twice} occur in more than one part")
            complete = False
        if [r["seed"] for r in rows] != list(SEEDS):
            why.append(f"V1: {key}: seeds {[r['seed'] for r in rows]}, not {list(SEEDS)}")
            complete = False
        for r in rows:
            for field in ("scores", "skills", "groups"):
                gone = sorted({f"{q}/{a}" for q in TRUTHS for a in ALL if a not in r.get(field, {}).get(q, {})})
                if gone:
                    why.append(f"V1: {key} seed {r['seed']}: missing from {field}: {gone[:6]}")
                    complete = False
            if (set(r.get("records", {})) != set(ARM_NAMES) or set(r.get("checks", {})) != set(LETTER.values())
                    or set(r.get("n_cells", {})) != set(ARM_NAMES)):
                why.append(f"V1: {key} seed {r['seed']}: records, checks or cell numbers missing")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"V1: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not every:
        return why or ["V1: no part"], dist
    if len({p["provenance"]["commit"] for p in every}) > 1:  # V2
        why.append("V2: the parts come from different commits")
    if not all(p["provenance"]["tree_clean"] for p in every):
        why.append("V2: a part comes from a dirty tree")
    if not commit or any(str(p["provenance"]["commit"]) != commit for p in every):
        why.append(f"V2: a part does not come from the run commit {commit!r} (the full hash)")
    if not complete:
        return why, dist
    for key in TASK_KEYS:
        ps, theirs = parts[key], h27_parts.get(key)
        for p in ps:
            why += table_checks(key, p, sha)
            if p.get("window_guard") != WINDOW_BYPASS:
                why.append(f"V6: {key}: the window guard's bypass is not recorded")
        if not theirs:
            why.append(f"V3: {key}: no committed row H27 part")
            continue
        ref = theirs[0]
        for p in ps:  # V3
            for q in TRUTHS:
                if not h27._same(p["truth_info"].get(q), ref["truth_info"].get(q), "cells"):
                    why.append(f"V3: {key}: truth {q} differs from row H27's (cells or core share)")
            for form in FORMS:
                if not h27._same(p["pred_info"].get(form), ref["pred_info"].get(form), "pool_cells"):
                    why.append(f"V3: {key}: {form}'s pool differs from row H27's (cells or core share)")
                f, g = p["forecast_targets"][form], ref["forecast_targets"][form]
                if not abs(f - g) <= h27.FACTOR_RTOL * abs(g):
                    why.append(f"V3: {key}: {form}'s forecast size differs from row H27's")
            if not p["v7"]["ok"]:
                why.append(f"V3: {key}: the reference cloud is not the cloud the size result measured")
            if len({json.dumps(q["table_sha256"]) for q in ps}) != 1:
                why.append(f"V6: {key}: the parts used different tables")
        rows = h26.merge_rows(ps)[0]
        trows = h26.merge_rows(theirs)[0]
        dist["V4"][key] = {}
        for form in FORMS:  # V4, the reproduction gate, per form
            w4, dist["V4"][key][LETTER[form]] = h28.gate(key, rows, trows, mine=name(form, CONTROL),
                                                         theirs=h27.arm_name(form, SIZING), truths=TRUTHS)
            why += [w.replace("V4: ", f"V4: {LETTER[form]} ", 1) for w in w4]
        why += identities(key, rows)
    return why, dist


# ------------------------------------------------------------------------------------------------ the tables
def series(rows: dict[str, list[dict]], keys, truths, nm: str, what: str = "score") -> np.ndarray:
    """Per seed: a score, skill or group skill (in board points) averaged over `truths`, then over `keys` seed by
    seed."""
    if what in GROUP_WEIGHT:
        def one(key, q):
            return np.array([100.0 * GROUP_WEIGHT[what] * r["groups"][q][nm][what] for r in rows[key]], float)
        return np.mean([np.mean([one(k, q) for q in truths], axis=0) for k in keys], axis=0)
    return h28.series(rows, keys, truths, nm, what)


def mean_series(rows, keys, truths, names, what: str = "score") -> np.ndarray:
    return np.mean([series(rows, keys, truths, n, what) for n in names], axis=0)


def compare(rows, keys, truths, a, b) -> dict:
    """a minus b (each a name or a tuple of names, averaged), paired by seed: score, every term, every group."""
    def s(x, what="score"):
        return mean_series(rows, keys, truths, (x,) if isinstance(x, str) else x, what)
    out = {"score": h28._stats(s(a) - s(b))}
    for m in (*TERMS, *SHAPE_TERMS):
        out[m] = float((s(a, m) - s(b, m)).mean())
    for g in GROUP_WEIGHT:
        out[f"points:{g}"] = float((s(a, g) - s(b, g)).mean())
    return out


UNITS = {**{k: (k,) for k in TASK_KEYS}, "fold E10.5": TASK_KEYS}


def pairs(form: str) -> list[tuple[str, object, object]]:
    n = lambda a: name(form, a)  # noqa: E731
    plc = tuple(n(p) for p in PLACEBOS)
    out = [(f"{n(a)} - {n(CONTROL)}", n(a), n(CONTROL)) for a in (CHORD0, *COUNTED)]
    out += [(f"{n(ARM_A2)} - {n(p)}", n(ARM_A2), n(p)) for p in PLACEBOS]
    out += [(f"{n(ARM_A2)} - mean placebo", n(ARM_A2), plc), (f"{n(ARM_A2)} - {n(ARM_A)}", n(ARM_A2), n(ARM_A))]
    return out


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: {lab: compare(rows, ks, ts, a, b) for f in FORMS for lab, a, b in pairs(f)}
                for u, ks in UNITS.items()} for w, ts in WORLDS}


def statements(tab: dict, form: str) -> dict:
    """S1 to S4 on one form, fixed in the row."""
    n = lambda a: name(form, a)  # noqa: E731

    def c(w, u, lab):
        return tab[w][u][lab]
    ws = [w for w, _ in WORLDS]
    a_ctl, a2_ctl = f"{n(ARM_A)} - {n(CONTROL)}", f"{n(ARM_A2)} - {n(CONTROL)}"
    s1 = all(c(w, "fold E10.5", a_ctl)["score"]["mean"] >= S1_FOLD
             and c(w, "fold E10.5", a_ctl)["score"]["ahead"] >= MIN_AHEAD
             and all(c(w, k, a_ctl)["score"]["mean"] >= S1_EACH for k in TASK_KEYS) for w in ws)
    s2 = all(c(w, "fold E10.5", f"{n(ARM_A2)} - {n(p)}")["score"]["mean"] >= S2_BAR
             and c(w, "fold E10.5", f"{n(ARM_A2)} - {n(p)}")["score"]["ahead"] >= MIN_AHEAD
             for w in ws for p in PLACEBOS)
    s3 = all(c(w, k, a2_ctl)["score"]["mean"] >= 0.0 for w in ws for k in TASK_KEYS)
    s4_parts = {}
    for w in ws:
        x = c(w, "fold E10.5", f"{n(ARM_A2)} - mean placebo")
        margin = x["score"]["mean"]
        share = sum(x[f"points:{g}"] for g in DE_GROUPS)
        s4_parts[w] = {"de": all(x[m] > 0 for m in DE_TERMS), "margin": margin, "de_groups_points": share,
                       "half": bool(margin > 0 and share >= 0.5 * margin)}
    s4 = all(v["de"] and v["half"] for v in s4_parts.values())
    return {"S1": bool(s1), "S2": bool(s2), "S3": bool(s3), "S4": bool(s4), "S4_parts": s4_parts,
            "all": bool(s1 and s2 and s3 and s4)}


def outcome(st: dict[str, dict], in_force: str, void: bool) -> str:
    """PASS: S1 to S4 on the form in force and on at least one other form. Otherwise, on the form in force: S1 not met
    NONE; S1 met and S2 not PLACEBO; S1 and S2 met and S4 not OFF_DE; anything else NONE."""
    if void or st is None:
        return "VOID"
    f = st[in_force]
    if f["all"] and any(st[o]["all"] for o in FORMS if o != in_force):
        return "PASS"
    if not f["S1"]:
        return "NONE"
    if not f["S2"]:
        return "PLACEBO"
    if not f["S4"]:
        return "OFF_DE"
    return "NONE"


# ------------------------------------------------------------------------------------------------ the form in force
FORM_WORDS = {"A": "full_iso_half", "B": "half_iso_half", "C": "core_iso_half", "full_iso_half": "full_iso_half",
              "half_iso_half": "half_iso_half", "core_iso_half": "core_iso_half"}
FORM_PHRASE = re.compile(r"named in this row before the read: ([A-Za-z_]+)")


def form_in_force(text: str | None) -> str | None:
    """The form the row names before the read (row H27's letter or name), or None."""
    if not text:
        return None
    m = FORM_PHRASE.findall(text)
    forms = {FORM_WORDS[w] for w in m if w in FORM_WORDS}
    return forms.pop() if len(forms) == 1 and len(m) == 1 else None


# ------------------------------------------------------------------------------------------------ the fixed words
FOLLOWS = {
    "VOID": ("VOID: {why}. No figure exists under row H30 and there is no rerun under it; the heart file in force "
             "stands."),
    "PASS": ("PASS: on the form in force ({form}) and on {others}, the chord at the practice tau is ahead of the "
             "control by the fixed bars, the test tau's chord beats each placebo, is not behind the control on either "
             "backtest, and gains on de_score and de_direction with at least half its margin over the placebos in "
             "expression_change and cell_state, in both worlds (one truth; an in-sample diagnostic). What follows: "
             "ONE priced question to the captain covering the ledger's widening of A15 to this outside input, the "
             "T2-heart amendments of X002 and X003, and the candidate beside the file in force with both worlds' "
             "practice figures as ranges and J6's caveat (tau ratios 2 and 1.33; the gain a range from zero to its "
             "practice value). The default is that the file in force stands; nothing is built before his word."),
    "PLACEBO": ("MATCHED BY ITS PLACEBO: on the form in force the chord is ahead of the control (S1) but does not beat "
                "each placebo by the bar (S2). H30 ends here: the captain is told in one line; nothing is asked."),
    "OFF_DE": ("OFF THE DE TERMS: on the form in force the chord is ahead of the control and of its placebos, but its "
               "margin over the placebos is not on de_score and de_direction, or not half in expression_change and "
               "cell_state (S4). H30 ends here: the captain is told in one line; nothing is asked."),
    "NONE": ("NONE: the fixed statements do not hold (S1 not met on the form in force, or S3 not met, or the pass not "
             "held on a second form). H30 ends here: the captain is told in one line; nothing is asked."),
}
L0 = ("Try H30 (run commit {commit}). A labelled in-sample diagnostic on one truth (E10.5, two backtests merged seed "
      "by seed); every figure is a practice reading on released stages and none is to be read beside a test score. "
      "Nothing is built and no file, rule, scope, size, count or depth changes. {bypass}.")


def follows(res: dict) -> str:
    o = res["outcome"]
    if o == "VOID":
        return FOLLOWS["VOID"].format(why="; ".join(res["why"]))
    if o == "PASS":
        f = res["form_in_force"]
        others = [LETTER[x] for x in FORMS if x != f and res["statements"][LETTER[x]]["all"]]
        return FOLLOWS["PASS"].format(form=LETTER[f], others="form " + " and form ".join(others))
    return FOLLOWS[o]


def reading(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], sha: dict, commit: str | None,
            in_force: str | None) -> dict:
    try:
        why, dist = void_reasons(parts, h27_parts, sha, commit)
    except (KeyError, TypeError, ValueError, IndexError) as e:  # a malformed part
        why, dist = [f"a part cannot be read ({type(e).__name__}: {e})"], {}
    if in_force not in FORMS:
        why = [*why, "the row names no form in force"]
    if why:
        res = {"void": True, "why": why, "guards": dist, "commit": commit, "outcome": "VOID",
               "form_in_force": in_force}
        res["follows"] = follows(res)
        return res
    rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
    tab = tables(rows)
    st = {f: statements(tab, f) for f in FORMS}
    res = {"void": False, "why": [], "row": ROW, "commit": commit, "form_in_force": in_force, "guards": dist,
           "tables": tab, "statements": {LETTER[f]: v for f, v in st.items()},
           "outcome": outcome(st, in_force, False),
           "chord_records": {k: parts[k][0]["table"] for k in TASK_KEYS},
           "window_guard": WINDOW_BYPASS}
    res["follows"] = follows(res)
    return res


def _st(s: dict) -> str:
    return f"{s['mean']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10 (se {s['se']:.3f})"


def text(res: dict) -> str:
    if res["void"]:
        return res["follows"] + "\n"
    lines = [L0.format(commit=res["commit"], bypass=WINDOW_BYPASS),
             f"Form in force (named in row {ROW} before the read): {LETTER[res['form_in_force']]}."]
    for key, g in res["guards"]["V4"].items():
        lines.append(f"GATE {key} (each form's control against row H27's committed form at f): " + "; ".join(
            f"{f} {v['mean']:.3f} / {v['seed']:.3f}" for f, v in g.items()) + f" (void beyond {GUARD_MEAN} / "
            f"{GUARD_SEED})")
    for f, s in res["statements"].items():
        lines.append(f"STATEMENTS form {f}: S1 {s['S1']}, S2 {s['S2']}, S3 {s['S3']}, S4 {s['S4']}"
                     + "".join(f"; {w}: DE {v['de']}, margin {v['margin']:+.3f}, in expression_change and cell_state "
                               f"{v['de_groups_points']:+.3f}" for w, v in s["S4_parts"].items()))
    lines += [f"OUTCOME: {res['outcome']}", res["follows"],
              "TABLES (decide nothing; paired by seed, mean [min, max] seeds ahead of 10 (se); per-term mean skill "
              "differences; per-group board points):"]
    for w, _ in WORLDS:
        for u in UNITS:
            for lab, c in res["tables"][w][u].items():
                lines.append(f"{TAG} {w} {u} {lab}: {_st(c['score'])}; "
                             + ", ".join(f"{m} {c[m]:+.4f}" for m in (*TERMS, *SHAPE_TERMS)) + "; "
                             + ", ".join(f"{g} {c[f'points:{g}']:+.3f}" for g in GROUP_WEIGHT))
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def cmd_gate(a) -> None:
    res = run_gate(Path(a.out), Path(a.data), Path(a.external))
    sys.exit(0 if res["outcome"] == "PASS" else 3)


def cmd_tables(a) -> None:
    run_tables(Path(a.out), Path(a.gate), Path(a.data), Path(a.external))


def cmd_run(a) -> None:
    refuse_run(Path(a.out) / PART)
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS}")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data), Path(a.tables))


def cmd_read(a) -> None:
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its read)")
    if GATE_ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: the gate's row {GATE_ROW} is not in experiments/harness/TRIES.md")
    out = Path(a.out)
    if (out / "reading.json").exists() or (out / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")
    if not gate_passed(COMMITTED_GATE):
        sys.exit(f"{ROW} is refused: the committed gate reading does not say PASS")
    if not committed_gate_ok(a.commit):
        sys.exit(f"{ROW} is refused: the committed gate reading was not made from a clean tree at an ancestor of the "
                 "run commit")
    in_force = form_in_force(row_text())
    if in_force is None:
        sys.exit(f"{ROW} is refused: row {ROW} does not name the form in force ('named in this row before the read: "
                 "<A, B or C>')")
    sha = json.loads(COMMITTED_SHA.read_text())
    res = reading(load_parts(Path(a.runs)), load_h27(Path(a.h27)), sha, a.commit, in_force)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "reading.txt").write_text(text(res))
    print(f"OUTCOME: {res['outcome']}")
    print(res["follows"])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--data", default=str(harness.DATA_DIR))
    g.add_argument("--external", default=str(EXTERNAL))
    t = sub.add_parser("tables")
    t.add_argument("--out", required=True)
    t.add_argument("--gate", required=True, help="the gate step's output folder (its gate.json must say PASS)")
    t.add_argument("--data", default=str(harness.DATA_DIR))
    t.add_argument("--external", default=str(EXTERNAL))
    r = sub.add_parser("run")
    r.add_argument("--task", required=True, help=f"one of {TASK_KEYS}")
    r.add_argument("--seeds", required=True, help="a range such as 0-1")
    r.add_argument("--out", required=True)
    r.add_argument("--tables", required=True, help="the tables step's output folder")
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--h27", default=str(h28.H27_PARTS))
    d.add_argument("--out", required=True)
    d.add_argument("--commit", required=True, help="the run commit (full hash) the parts must come from")
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "tables": cmd_tables, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
