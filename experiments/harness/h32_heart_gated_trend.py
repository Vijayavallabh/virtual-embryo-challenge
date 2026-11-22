#!/usr/bin/env python
"""Try H32 of experiments/harness/TRIES.md: the heart's ladder trend (row H28's trend_s1) in cardiomyocytes and in
endothelial cells (Endocardium plus Endothelium, b0_type 'Endothelial') gated gene by gene by the sign of the outside
late chord (row H30's tables), on form C (the full trim), the file in force. A score-blind round's specification (the
cold challenger's final entry). An expression step only: no cell, position, size, scope, count, depth or composition
changes. The outside data only choose genes; they set no size.

The trend. Exactly as h28.trend_expression(ctx, 1) builds it: harness.b0_trend_factors on the context's draw, typed by
obs['merged'], damping harness.B0_DAMPING, seed [ctx.seed, 3]. For k in ('CM', 'Endothelial') the base factor b_k is the
mapping's logf[k] ("mapping" case) or, for a lineage the mapping leaves out (usable at start and the reference only),
B0_DAMPING x (y_ref[k] - y_start[k]) with y from h28.lineage_log_means ("fallback" case), as h28's placebo builds its
`full` dict. A lineage momentum_shift skips for too few cells stays skipped ("skipped"); one the draw does not hold is
"absent". The case is recorded per lineage and seed.

The chord. h30.build_table(key, ...) unchanged (the corrected reader of rows N56 and H30), one table per backtest; the
test table is the one-day backtest's table (its anchors are the test's: JAX 22 and Feng E9.5 early, JAX 54-56 and Feng
E14.5 late, classifier on MERFISH E9.5 labels). Pairing: 'CM' with chord lineage 'CM', 'Endothelial' with 'Endo'. A
paired chord lineage absent from a table leaves that type ungated in every arm (recorded).

The gate. T_D = 0.25. E_k = {g : b_k[g] != 0 and |D_k[g]| >= T_D}; M_k[g] = 1 on E_k where sign(D) == sign(b), else 0.
G's mapping is the ungated mapping with mapping[k] = b_k * M_k for the paired k; every other key unchanged; it enters
at logf_by_type in the momentum_shift call.

Arms (per backtest, per seed; everything else is row H28's: draw, coordinates, row H27's f size step, row H25's half,
the checked builder):
  control  row H27's C at f (reference)
  U        the ungated trend through this script's own path with every mask one (reference; V4 checks it equals
           h28.trend_expression(ctx, 1) bit for bit)
  G        the gate
  R1-R3    per k, E_k split into ten equal-count bins of the lineage's reference-stage log(mean CP10k + 0.5)
           (h30.reference_basis, whole reference stage); in each bin as many genes as G keeps there, drawn uniformly
           without replacement, seed 20261130 + j + crc32(k); kept genes get b_k, the others zero
  F        the sign-flipped gate: E_k genes whose sign disagrees keep b_k; others zero
  P1-P3    G's gated vector per k deranged within ten equal-count bins of the same level over all genes
           (h28.derangement, seed 20261126 + j); no gene keeps its own factor

Steps (one per call, run_h32.sh): tests; tables (row H30's tables, or H30's own committed files where its tables step
wrote them at the same reader); gate (G0, target-free, not a try: the agreement a at the test pair against a null of
1000 binned permutations); run (one backtest, a seed range, one JSON); read (the void conditions, the fixed statements
S1 to S6, AHEAD or NOT AHEAD; prints the outcome and what follows only). The harness's external-window guard is
bypassed by name, as row H30 does, and every part records it.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
import zlib
from dataclasses import replace
from pathlib import Path

import numpy as np

from vec import harness, prep, provenance, t1_scope, t2_build
from vec.adoption import SEEDS
from vec.boards import load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parents[2]  # this checkout


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h30 = _load("h30_heart_late_chord", HERE / "experiments" / "harness" / "h30_heart_late_chord.py")
h28 = h30.h28  # the one h28 module h30 loaded (trend_expression, derangement, lineage_log_means, gate, series)
h27, h26, h25, payoff = h28.h27, h28.h26, h28.h25, h28.payoff

ROW = "H32"
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
REGISTER = HERE / "EXTERNAL_DATA.md"
REGISTER_ROWS = ("X002", "X003")
AMENDMENT = "T2-heart amendment for row H32"  # must open a sentence or a cell of each register row's line
AXIS = "h32-heart-gated-trend"
PART = "h32.json"
SOURCE = "experiments/harness/h32_heart_gated_trend.py (try H32)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h32_heart_gated_trend"
COMMITTED_GATE = RESULTS / "gate" / "gate.json"
COMMITTED_SHA = RESULTS / "tables_sha256.json"
H28_PARTS = HERE / "experiments" / "harness" / "results" / "h28_heart_trend_strength" / "parts"
# the reader code the tables depend on (a table of row H30 is reused only where each blob is the same at its commit)
READER_FILES = ("experiments/harness/h30_heart_late_chord.py", "experiments/external/feng2022_geo.py",
                "experiments/external/feng_relabel.py", "experiments/external/filterio.py",
                "experiments/external/windows.py")

# ------------------------------------------------------------------------------------------------ fixed numbers
T_D = 0.25  # the scorer's |lfc| threshold for a DE gene
PAIR = {"CM": "CM", "Endothelial": "Endo"}  # b0_type -> chord lineage
N_BINS = 10
R_SEED, P_SEED, G0_SEED = 20261130, 20261126, 20261122
N_R = N_P = 3
G0_PERM, G0_MIN_GENES, G0_MIN_A = 1000, 50, 0.60
TASK_KEYS = h30.TASK_KEYS  # T2h/<=E9.5->E10.5, T2h/<=E8.75->E10.5
TEST_TABLE_KEY = "T2h/<=E9.5->E10.5"  # the test uses the one-day backtest's anchors and classifier
FORM = "core_iso_half"  # form C
SIZING = "f"
H27_ARM = h27.arm_name(FORM, SIZING)
CONTROL, U, G, F = "control", "U", "G", "F"
RS = tuple(f"R{j}" for j in range(1, N_R + 1))
PS = tuple(f"P{j}" for j in range(1, N_P + 1))
COUNTED = (G, *RS, F, *PS)
BUILT = (CONTROL, U, *COUNTED)
ALL = (harness.FLOOR, *BUILT)
OWN, CLEAN = h26.OWN, h26.CLEAN
TRUTHS = (OWN, *CLEAN)
WORLDS = (("own", (OWN,)), ("CLEAN", CLEAN))
GUARD_MEAN, GUARD_SEED = h25.GUARD_MEAN, h25.GUARD_SEED  # 0.1 and 0.5 board points
U_SCORE_TOL = 0.01  # U against row H28's committed trend_s1, per seed and truth
TOL = 1e-9
SHAPE_TERMS = h28.SHAPE_TERMS
DE_TERMS = ("de_score", "de_direction")
TERMS = h30.TERMS
DE_GROUPS = ("expression_change", "cell_state")
BAR, MIN_AHEAD, F_MIN_AHEAD = 0.3, 9, 8
CASES = ("mapping", "fallback", "skipped", "absent")
TEST_ROLES = {"earliest": "T2/E8.25_late.h5ad", "earlier2": "T2/E8.5.h5ad", "earlier": "T2/E8.75.h5ad",
              "start": "T2/E9.5.h5ad", "ref": "T2/E10.5.h5ad"}  # row N27's inputs
TEST_BOARD = "T2:heart:test_extrap"
PRACTICE_BUDGET = h28.PRACTICE_BUDGET
TEST_BUDGET = h28.TEST_BUDGET
TAG = "[practice diagnostic, released stages, one truth, in-sample]"
WINDOW_BYPASS = h30.WINDOW_BYPASS


# ------------------------------------------------------------------------------------------------ refusals
def tries_rows(path: Path = TRIES) -> set[str]:
    return h27.tries_rows(path)


def register_missing(path: Path = REGISTER, ids=REGISTER_ROWS, words: str = AMENDMENT) -> list[str]:
    """The register rows among `ids` whose line does not carry `words` opening a sentence or a cell."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    out = []
    for x in ids:
        ln = next((ln for ln in lines if ln.startswith(f"| {x} ")), None)
        if ln is None or not h30._opens_sentence(ln, words):
            out.append(x)
    return out


def refuse_rows_and_register(tries: Path = TRIES, register: Path = REGISTER) -> None:
    if ROW not in tries_rows(tries):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    missing = register_missing(register)
    if missing:
        sys.exit(f"{ROW} is refused: EXTERNAL_DATA.md rows {missing} do not carry '{AMENDMENT}' (it comes before any "
                 "outside file is read)")


def _clean_or_exit() -> str:
    commit, unclean = provenance.git_state(HERE)
    if unclean is None or unclean:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    return commit


def refuse_common(out_files) -> str:
    refuse_rows_and_register()
    commit = _clean_or_exit()
    for p in out_files:
        if Path(p).exists():
            sys.exit(f"{ROW} is refused: {p} exists")
    return commit


def gate_passed(path: Path, commit: str | None = None) -> bool:
    """G0's reading at `path` says PASS (and, given a commit, was made at it from a clean tree)."""
    try:
        g = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return False
    if g.get("outcome") != "PASS" or g.get("axis") != AXIS + ":gate":
        return False
    if commit is not None:
        return g.get("provenance", {}).get("commit") == commit and g.get("provenance", {}).get("tree_clean") is True
    return True


def committed_gate_ok(commit: str, path: Path = COMMITTED_GATE) -> bool:
    """The committed G0 reading says PASS and was made from a clean tree at a strict ancestor of `commit`."""
    if not gate_passed(path):
        return False
    prov = json.loads(Path(path).read_text()).get("provenance", {})
    return prov.get("tree_clean") is True and h30.is_ancestor(str(prov.get("commit", "")), str(commit), HERE)


# ------------------------------------------------------------------------------------------------ the gate's sets
def chord_vector(table: dict, genes, lineage: str) -> np.ndarray | None:
    """The table's D of `lineage` in `genes` order (0 for a gene the table does not hold), or None when the table
    does not move that lineage."""
    if lineage not in table["lineages"]:
        return None
    pos = {g: i for i, g in enumerate(table["genes"])}
    d = np.asarray(table["lineages"][lineage], dtype=np.float64)
    idx = np.array([pos.get(g, -1) for g in genes], dtype=np.int64)
    return np.where(idx >= 0, d[np.maximum(idx, 0)], 0.0)


def gate_sets(b: np.ndarray, d: np.ndarray) -> dict[str, np.ndarray]:
    """E (eligible: b != 0 and |D| >= T_D), M (E and the signs agree) and Fl (E and the signs disagree)."""
    b, d = np.asarray(b, np.float64), np.asarray(d, np.float64)
    e = (b != 0) & (np.abs(d) >= T_D)
    m = e & (np.sign(d) == np.sign(b))
    return {"E": e, "M": m, "Fl": e & ~m}


def level_bins(idx: np.ndarray, basis: np.ndarray, bins: int = N_BINS) -> list[np.ndarray]:
    """`idx` split into `bins` equal-count bins by `basis` (rank, stable)."""
    idx = np.asarray(idx, dtype=np.int64)
    return np.array_split(idx[np.argsort(np.asarray(basis, np.float64)[idx], kind="stable")], bins)


def random_gate(e: np.ndarray, m: np.ndarray, basis: np.ndarray, k: str, j: int) -> np.ndarray:
    """R_j's kept genes for type k: per level bin of E, as many genes as M keeps there, uniformly without replacement
    from that bin, seed R_SEED + j + crc32(k)."""
    rng = np.random.default_rng(R_SEED + j + zlib.crc32(k.encode()))
    out = np.zeros(len(e), dtype=bool)
    for b in level_bins(np.flatnonzero(e), basis):
        n = int(np.asarray(m)[b].sum())
        if n:
            out[rng.choice(b, size=n, replace=False)] = True
    return out


def bin_counts(mask: np.ndarray, e: np.ndarray, basis: np.ndarray) -> list[int]:
    return [int(np.asarray(mask)[b].sum()) for b in level_bins(np.flatnonzero(e), basis)]


def agreement(b: np.ndarray, d: np.ndarray) -> tuple[int, float]:
    s = gate_sets(b, d)
    n = int(s["E"].sum())
    return n, (float(s["M"].sum() / n) if n else float("nan"))


def g0_lineage(b: np.ndarray, d: np.ndarray, basis: np.ndarray, lineage: str, n_perm: int = G0_PERM) -> dict:
    """G0 for one lineage: n = |E|, a = the share of E whose signs agree, and a's null: D permuted among E's genes
    within ten equal-count bins of `basis` (E fixed), n_perm times, seed G0_SEED + crc32(lineage)."""
    s = gate_sets(b, d)
    e = np.flatnonzero(s["E"])
    n, a = agreement(b, d)
    rec = {"lineage": lineage, "n": n, "a": a, "kept": int(s["M"].sum())}
    if n == 0:
        return {**rec, "null_p95": float("nan"), "null_mean": float("nan")}
    rng = np.random.default_rng(G0_SEED + zlib.crc32(lineage.encode()))
    sb, sd = np.sign(np.asarray(b, float)[e]), np.asarray(d, float)[e]
    pos = {g: i for i, g in enumerate(e)}
    bins = [np.array([pos[g] for g in bb], dtype=np.int64) for bb in level_bins(e, basis)]
    null = np.empty(n_perm)
    for i in range(n_perm):
        p = np.arange(n)
        for bb in bins:
            p[bb] = rng.permutation(bb)
        null[i] = float(np.mean(np.sign(sd[p]) == sb))
    return {**rec, "null_p95": float(np.percentile(null, 95)), "null_mean": float(null.mean())}


def g0_outcome(cm: dict) -> str:
    ok = cm.get("n", 0) >= G0_MIN_GENES and cm.get("a", float("nan")) >= G0_MIN_A and cm["a"] > cm["null_p95"]
    return "PASS" if ok else "FAIL"


# ------------------------------------------------------------------------------------------------ the trend
def trend_parts(ctx: harness.BuildContext) -> dict:
    """h28.trend_expression(ctx, 1)'s inputs, line by line: the draw, the typed cells, the mapping and its record."""
    harness._b0_trend_task_check(ctx.task)
    pool, rows, _, coords = harness.draw_copy(ctx)
    like = h28._typed(pool, "end")[rows]
    genes = list(like.var_names)
    types = sorted(set(np.asarray(like.obs["b0_type"]).astype(str)))
    logf, rec = harness.b0_trend_factors(ctx, genes, types)
    return {"pool": pool, "rows": rows, "coords": coords, "start": ctx.inputs["start"], "like": like, "genes": genes,
            "types": types, "logf": logf, "rec": rec}


def base_factors(ctx: harness.BuildContext, parts: dict, types=tuple(PAIR)) -> tuple[dict, dict]:
    """b_k for each paired type and its case: the mapping's factor, or B0's own (fallback) as h28's placebo builds
    it, or none (skipped for too few cells; absent from the draw)."""
    logf, rec, genes = parts["logf"], parts["rec"], parts["genes"]
    fallback = {r["lineage"] for r in rec["lineages"] if r["own_factor"]}
    b, cases = {}, {}
    for k in types:
        if k in logf:
            b[k], cases[k] = np.asarray(logf[k], dtype=np.float64), "mapping"
        elif k in fallback:
            y_ref = h28.lineage_log_means(ctx, ctx.task.reference_role, genes, [k])
            y_start = h28.lineage_log_means(ctx, "start", genes, [k])
            b[k], cases[k] = 1.0 * (harness.B0_DAMPING * (y_ref[k] - y_start[k])), "fallback"
        else:
            cases[k] = "skipped" if k in parts["types"] else "absent"
    return b, cases


def express(ctx: harness.BuildContext, parts: dict, replace_by_type: dict[str, np.ndarray]):
    """The trend's expression with the mapping's entries for the named types replaced (a fallback type named here
    takes the given factor in place of B0's own); returns the stage through harness.finish and the shift report."""
    mapping = {k: 1.0 * v for k, v in parts["logf"].items()}
    mapping.update({k: np.asarray(v, dtype=np.float64) for k, v in replace_by_type.items()})
    shifted, report = t1_scope.momentum_shift(h28._typed(parts["start"], "start"), h28._typed(parts["pool"], "end"),
                                              parts["like"], harness.B0_DAMPING * 1.0, (), key="b0_type",
                                              seed=[ctx.seed, 3], logf_by_type=mapping)
    st = harness.finish(ctx, parts["pool"], dense(shifted[:, ctx.board.genes].X), parts["coords"])
    return st, report


def arm_vectors(b: dict, chords: dict, basis: dict) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    """Every arm's replacement vectors per gated type (b_k times the arm's mask, or P's derangement of G's vector),
    and the record of the sets: U ones; G M; R_j the random gate; F the flipped set; P_j G deranged. A type with no b
    or no chord is in no arm's replacement (it keeps its ungated factor, or stays skipped)."""
    out = {a: {} for a in (U, G, *RS, F, *PS)}
    rec = {}
    for k, bk in b.items():
        d = chords.get(k)
        if d is None:
            rec[k] = {"gated": False}
            continue
        s = gate_sets(bk, d)
        lvl = basis[k]
        out[U][k] = bk * np.ones(len(bk))
        g = bk * s["M"]
        out[G][k] = g
        r = {}
        for j, name in enumerate(RS, start=1):
            m = random_gate(s["E"], s["M"], lvl, k, j)
            out[name][k] = bk * m
            r[name] = {"inside_E": bool(not np.any(m & ~s["E"])),
                       "bin_counts_match": bin_counts(m, s["E"], lvl) == bin_counts(s["M"], s["E"], lvl),
                       "kept": int(m.sum())}
        out[F][k] = bk * s["Fl"]
        fixed = {}
        for j, name in enumerate(PS, start=1):
            p = h28.derangement(lvl, k, bins=N_BINS, seed=P_SEED + j)
            out[name][k] = g[p]
            fixed[name] = int((p == np.arange(len(p))).sum())
        rec[k] = {"gated": True, "eligible": int(s["E"].sum()), "kept": int(s["M"].sum()), "flipped": int(s["Fl"].sum()),
                  "f_complement": bool(np.array_equal(s["Fl"] | s["M"], s["E"]) and not np.any(s["Fl"] & s["M"])),
                  "R": r, "P_fixed_points": fixed, "kept_bin_counts": bin_counts(s["M"], s["E"], lvl)}
    return out, rec


def type_basis(whole_ref, genes, types=tuple(PAIR)) -> dict[str, np.ndarray]:
    """Each paired type's level: its chord lineage's log(mean CP10k + 0.5) in the whole reference stage (released
    labels; h30.reference_basis)."""
    lb = h30.reference_basis(whole_ref, list(genes), sorted({PAIR[k] for k in types}))
    return {k: lb[PAIR[k]] for k in types}


def build_arms(ctx: harness.BuildContext, size: dict, table: dict, whole_ref) -> tuple[dict, dict]:
    """Every arm of one seed on form C's context, each then sized by row H27's f step, halved and checked; and the
    checks of V4 to V6."""
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    parts = trend_parts(ctx)
    b, cases = base_factors(ctx, parts)
    chords = {k: chord_vector(table, parts["genes"], PAIR[k]) for k in b}
    basis = type_basis(whole_ref, parts["genes"])
    vecs, sets = arm_vectors(b, chords, basis)
    pre, reports = {CONTROL: harness._copy(ctx)}, {}
    for a in (U, *COUNTED):
        pre[a], rep = express(ctx, parts, vecs[a])
        reports[a] = sorted(r["type"] for r in rep.rows if r["status"] == "shifted")
    ref_u, _, ref_rows, _ = h28.trend_expression(ctx, 1.0)  # row H28's trend_s1 on the same context (V4)
    ctrl = pre[CONTROL]
    checks = {"u_equals_h28_trend_s1": bool(np.array_equal(pre[U].X, ref_u.X)
                                            and np.array_equal(pre[U].coords, ref_u.coords)),
              "same_draw": bool(np.array_equal(ref_rows, parts["rows"])),
              "pre_size_coords_equal": bool(all(np.array_equal(v.coords, ctrl.coords) for v in pre.values())),
              "cases": cases, "sets": sets,
              "chord_lineages": {k: (PAIR[k] if chords.get(k) is not None else None) for k in b}}
    out, recs = {}, {}
    for name, st in pre.items():
        X, C, srec = h27.size_arm(st.X, st.coords, FORM, SIZING, {}, size, {}, seed)
        sub, rep = t2_build.build_submission(ctx.board, X, C, genes=list(ctx.board.genes))
        if not rep.ok:
            raise ValueError(f"{ctx.task.key} {name}: the built file fails its check: {rep.errors}")
        out[name] = Stage.from_anndata(sub, ctx.board.genes, need_coords=True)
        recs[name] = {"size": srec, "shifted": reports.get(name)}
    checks["post_size_coords_equal"] = bool(all(np.array_equal(v.coords, out[CONTROL].coords) for v in out.values()))
    checks["n_cells"] = {k: int(v.n) for k, v in out.items()}
    checks["trend"] = {"budget_days": float(parts["rec"]["budget_days"]), "stages": parts["rec"]["stages"],
                       "eps": float(parts["rec"]["eps"]), "n_trend": int(parts["rec"]["n_trend"])}
    return out, {"records": recs, "checks": checks}


# ------------------------------------------------------------------------------------------------ the tables step
def reader_same(commit: str, repo: Path = HERE) -> bool:
    """Each reader file's blob at `commit` is its blob at HEAD."""
    for f in READER_FILES:
        a = subprocess.run(["git", "rev-parse", f"{commit}:{f}"], cwd=repo, capture_output=True, text=True)
        b = subprocess.run(["git", "rev-parse", f"HEAD:{f}"], cwd=repo, capture_output=True, text=True)
        if a.returncode or b.returncode or a.stdout.strip() != b.stdout.strip():
            return False
    return True


def h30_reusable(h30_tables: Path | None, committed: Path = h30.COMMITTED_SHA, same=reader_same) -> dict | None:
    """Row H30's own tables, where its tables step wrote them and their sha256 is committed: each file's sha256 is the
    committed one and each was made at a commit whose reader is this commit's. Returns {name: bytes} or None."""
    if h30_tables is None or not Path(committed).exists():
        return None
    try:
        want = json.loads(Path(committed).read_text())
    except (OSError, ValueError):
        return None
    got = {}
    for key in TASK_KEYS:
        nm = f"{h30.TRANSITIONS[key]['short']}.json"
        p = Path(h30_tables) / nm
        if not p.exists():
            return None
        by = p.read_bytes()
        if want.get(nm) != h30.sha256_bytes(by):
            return None
        try:
            commit = json.loads(by)["provenance"]["commit"]
        except (ValueError, KeyError, TypeError):
            return None
        if not same(str(commit)):
            return None
        got[nm] = by
    return got


def run_tables(out: Path, data_dir: Path, external: Path, h30_tables: Path | None) -> dict:
    out = Path(out)
    files = [out / f"{h30.TRANSITIONS[k]['short']}.json" for k in TASK_KEYS] + [out / "tables_sha256.json",
                                                                                 out / "tables_origin.json"]
    commit = refuse_common(files)
    out.mkdir(parents=True, exist_ok=True)
    reused = h30_reusable(h30_tables)
    sums = {}
    if reused is not None:
        for nm, by in reused.items():
            (out / nm).write_bytes(by)
            sums[nm] = h30.sha256_bytes(by)
        origin = {"reused_from_row_H30": True, "from": str(h30_tables)}
    else:
        h30.need_rdata()
        panel = list(load_boards()[payoff.heart_task(TASK_KEYS[0]).panel_board].genes)
        cache: dict = {}
        for key in TASK_KEYS:
            by = h30.table_bytes(h30.build_table(key, data_dir, external, panel, cache))
            nm = f"{h30.TRANSITIONS[key]['short']}.json"
            (out / nm).write_bytes(by)
            sums[nm] = h30.sha256_bytes(by)
            print(f"tables {h30.TRANSITIONS[key]['short']} written", flush=True)
        origin = {"reused_from_row_H30": False}
    origin.update(commit=commit, test_table=f"{h30.TRANSITIONS[TEST_TABLE_KEY]['short']}.json")
    (out / "tables_origin.json").write_text(json.dumps(origin, indent=1, sort_keys=True) + "\n")
    (out / "tables_sha256.json").write_text(json.dumps(sums, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out / 'tables_sha256.json'}", flush=True)
    return sums


def load_table(path: Path, committed: Path) -> tuple[dict, str]:
    return h30.load_table(Path(path), Path(committed))


# ------------------------------------------------------------------------------------------------ G0
def _test_view(board_key: str) -> harness.TaskView:
    days = {r: float(harness.STAGE_DAYS[f]) for r, f in TEST_ROLES.items()}
    return harness.TaskView(key="test:T2h", task="T2", kind="extrap", setting="heart", test_board=board_key,
                            roles=tuple(TEST_ROLES), role_days=days,
                            role_stems={r: harness.stem(f) for r, f in TEST_ROLES.items()}, reference_role="ref",
                            target_day=float("nan"))  # never read by the trend factor


def pair_values(ctx: harness.BuildContext, pool, whole_ref, table: dict) -> dict:
    """G0's statistic on one pair: b from the pool's own types (no draw), D from the table, per paired type."""
    genes = list(pool.var_names)
    types = sorted(set(np.asarray(pool.obs[harness.T2_LINEAGE_KEY]).astype(str)))
    logf, rec = harness.b0_trend_factors(ctx, genes, types)
    parts = {"logf": logf, "rec": rec, "genes": genes, "types": types}
    b, cases = base_factors(ctx, parts)
    basis = type_basis(whole_ref, genes)
    out = {}
    for k in PAIR:
        d = chord_vector(table, genes, PAIR[k])
        if k not in b or d is None:
            out[PAIR[k]] = {"lineage": PAIR[k], "case": cases[k], "chord": d is not None, "n": 0, "a": float("nan"),
                            "null_p95": float("nan"), "null_mean": float("nan")}
            continue
        out[PAIR[k]] = {**g0_lineage(b[k], d, basis[k], PAIR[k]), "case": cases[k], "chord": True}
    return out


def _pool_c(whole):
    pcore = np.asarray(whole.obs["core"], dtype=bool)
    PC = np.asarray(whole.obsm[COORD_KEY], dtype=float)[:, :3]
    return whole[h26.pool_masks(PC, pcore)[FORM]].copy()


def n56_outcome(path: Path = h30.COMMITTED_GATE) -> str:
    """Row N56's committed outcome word, reported beside G0; it decides nothing."""
    try:
        return str(json.loads(Path(path).read_text()).get("outcome", "not recorded"))
    except (OSError, ValueError):
        return "not committed"


def run_gate(out: Path, tables: Path, data_dir: Path) -> dict:
    out = Path(out)
    commit = refuse_common([out / "gate.json", out / "gate.txt"])
    sha = Path(tables) / "tables_sha256.json"
    test_table, test_sha = load_table(Path(tables) / f"{h30.TRANSITIONS[TEST_TABLE_KEY]['short']}.json", sha)
    if test_table.get("task") != TEST_TABLE_KEY:
        sys.exit(f"{ROW} is refused: the test table is not the one-day backtest's")
    try:
        origin = json.loads((Path(tables) / "tables_origin.json").read_text())
    except (OSError, ValueError):
        sys.exit(f"{ROW} is refused: the tables step's origin record is missing")
    boards = load_boards()
    board = boards[TEST_BOARD]
    view = _test_view(board.key)
    inputs = {r: load_stage(Path(data_dir) / f, genes=board.genes, drop_duplicates=True) for r, f in TEST_ROLES.items()}
    whole = inputs["ref"]
    pool = _pool_c(whole)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs={**inputs, "ref": pool},
                               seed=harness.BUILD_SEED_OFFSET)
    test = pair_values(ctx, pool, whole, test_table)
    outcome = g0_outcome(test["CM"])
    backtests = {}
    for key in TASK_KEYS:  # reported, deciding nothing
        tb, _ = load_table(Path(tables) / f"{h30.TRANSITIONS[key]['short']}.json", sha)
        t = payoff.heart_task(key)
        genes = boards[t.panel_board].genes
        ins = {r: load_stage(Path(data_dir) / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
        p = _pool_c(ins["ref"])
        bctx = harness.BuildContext(task=harness.TaskView.of(t), board=boards[t.panel_board],
                                    inputs={**ins, "ref": p}, seed=harness.BUILD_SEED_OFFSET)
        backtests[key] = pair_values(bctx, p, ins["ref"], tb)
    res = {"axis": AXIS + ":gate", "row": ROW, "outcome": outcome, "test_pair": test, "backtests": backtests,
           "test_table_sha256": test_sha, "tables_origin": origin, "n56_outcome_reported": n56_outcome(),
           "fixed": {"t_d": T_D, "min_genes": G0_MIN_GENES, "min_a": G0_MIN_A, "perm": G0_PERM, "bins": N_BINS,
                     "seed": G0_SEED},
           "note": "target-free: compares no predictor with a fold target; counts no try",
           "source": SOURCE, "provenance": {**payoff.provenance(), "commit_at_start": commit}}
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "gate.txt").write_text(gate_text(res))
    print(f"gate {outcome}", flush=True)
    return res


def _g0_line(where: str, r: dict) -> str:
    return (f"{where} {r['lineage']}: case {r.get('case')}; |E| {r['n']}; a {r['a']:.3f}; null 95th percentile "
            f"{r['null_p95']:.3f} (mean {r['null_mean']:.3f})")


def gate_text(res: dict) -> str:
    lines = [f"Try H32, G0 (commit {res['provenance']['commit']}): target-free, counts no try. OUTCOME: "
             f"{res['outcome']}. PASS needs, at the test pair for CM, |E| at least {G0_MIN_GENES}, a at least "
             f"{G0_MIN_A} and a strictly above its null's 95th percentile; otherwise nothing further runs or counts."]
    lines += [_g0_line("test pair", r) for r in res["test_pair"].values()]
    lines += [_g0_line(f"backtest {k} (reported, decides nothing)", r) for k, v in res["backtests"].items()
              for r in v.values()]
    lines.append(f"Row N56's committed outcome, reported beside G0 (decides nothing): {res['n56_outcome_reported']}.")
    lines.append(f"Tables reused from row H30: {res['tables_origin'].get('reused_from_row_H30')}.")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ the run
def refuse_run(path: Path) -> str:
    commit = refuse_common([path])
    if not gate_passed(COMMITTED_GATE):
        sys.exit(f"{ROW} is refused: the committed G0 reading does not say PASS (nothing further runs)")
    if not committed_gate_ok(commit):
        sys.exit(f"{ROW} is refused: the committed G0 reading was not made from a clean tree at an ancestor of this "
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
    table, table_sha = load_table(Path(tables) / f"{h30.TRANSITIONS[key]['short']}.json", COMMITTED_SHA)
    if table.get("task") != key:
        sys.exit(f"{ROW} is refused: the table is not {key}'s")
    t = payoff.heart_task(key)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    geoms, _ = h27.input_geometries(t, data_dir)
    size = h27.analogue_size_result(geoms, float(harness.STAGE_DAYS[t.target]))
    ftargets = h27.forecast_targets(size)
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    tm = h26.truth_masks(TC, tcore)
    truths = {q: Stage.from_anndata(ta[tm[q]], genes, need_coords=True) for q in TRUTHS}
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    whole = inputs["ref"]
    pool = _pool_c(whole)
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        ctx = harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pool},
                                   seed=harness.BUILD_SEED_OFFSET + s)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            preds, extra = build_arms(ctx, size, table, whole)
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": extra["checks"]["n_cells"],
               "checks": extra["checks"], "records": extra["records"], "scores": {}, "skills": {}, "groups": {}}
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
    res = {"label": "try H32: the heart's ladder trend gated by the sign of the outside late chord, form C; a "
                    "labelled in-sample diagnostic on one truth, not an adoption result", "axis": AXIS, "task": t.key,
           "arms": list(BUILT), "truth_values": list(TRUTHS), "forecast_targets": ftargets,
           "truth_seed": payoff.TRUTH_SEED, "half_seed_base": h25.HALF_SEED_BASE, "table_sha256": table_sha,
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


# ------------------------------------------------------------------------------------------------ the void
def u_scores(key: str, rows: list[dict], h28_rows: list[dict]) -> list[str]:
    """V4's second half: U's board score within U_SCORE_TOL of row H28's committed trend_s1, per seed and truth."""
    other = {r["seed"]: r for r in h28_rows}
    why = []
    for r in rows:
        for q in TRUTHS:
            try:
                d = abs(r["scores"][q][U] - other[r["seed"]]["scores"][q][h28.arm(1.0)])
            except KeyError:
                why.append(f"V4: {key} seed {r['seed']} {q}: no committed row H28 trend_s1 score to compare")
                continue
            if not d <= U_SCORE_TOL:
                why.append(f"V4: {key} seed {r['seed']} {q}: U is {d:.4f} from row H28's committed trend_s1")
    return why


def identities(key: str, rows: list[dict]) -> list[str]:
    """V4 (U's X is h28's trend_s1 bit for bit), V5 (draw, coordinates, cell number, shape terms) and V6's per-seed
    parts (R inside E with G's bin counts, F the complement, P no fixed point, every paired type's case recorded)."""
    why = []
    for r in rows:
        c = r["checks"]
        tag = f"{key} seed {r['seed']}"
        if c.get("u_equals_h28_trend_s1") is not True:
            why.append(f"V4: {tag}: U's X is not h28.trend_expression(ctx, 1)'s bit for bit")
        for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
            if c.get(k) is not True:
                why.append(f"V5: {tag}: {k.replace('_', ' ')} fails")
        if set(r["n_cells"]) != set(BUILT) or len(set(r["n_cells"].values())) != 1:
            why.append(f"V5: {tag}: the arms' cell numbers differ")
        for q in TRUTHS:
            for a in BUILT:
                for m in SHAPE_TERMS:
                    if abs(r["skills"][q][a][m] - r["skills"][q][CONTROL][m]) > TOL:
                        why.append(f"V5: {tag} {q}: {a}'s {m} is not the control's")
        cases = c.get("cases", {})
        if set(cases) != set(PAIR) or any(v not in CASES for v in cases.values()):
            why.append(f"V6: {tag}: a paired type's case (mapping, fallback, skipped, absent) is not recorded")
        for k, s in c.get("sets", {}).items():
            if not s.get("gated"):
                continue
            for nm in RS:
                if s["R"][nm]["inside_E"] is not True:
                    why.append(f"V6: {tag} {k}: {nm} draws outside E")
                if s["R"][nm]["bin_counts_match"] is not True:
                    why.append(f"V6: {tag} {k}: {nm} misses G's per-bin count")
            if s.get("f_complement") is not True:
                why.append(f"V6: {tag} {k}: F is not G's complement in E")
            if any(v for v in s.get("P_fixed_points", {}).values()):
                why.append(f"V6: {tag} {k}: a P arm leaves a gene on its own factor")
    return why


def void_reasons(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], h28_parts: dict[str, list[dict]],
                 sha: dict, commit: str | None, gate_ok: bool) -> tuple[list[str], dict]:
    why: list[str] = []
    dist: dict = {"V3": {}}
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
            if set(r.get("records", {})) != set(BUILT) or "checks" not in r or set(r.get("n_cells", {})) != set(BUILT):
                why.append(f"V1: {key} seed {r['seed']}: records, checks or cell numbers missing")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"V1: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not gate_ok:
        why.append("V6: G0 is not committed as PASS")
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
        ps = parts[key]
        for p in ps:
            why += h30.table_checks(key, p, sha)
            if p.get("window_guard") != WINDOW_BYPASS:
                why.append(f"V6: {key}: the window guard's bypass is not recorded")
        if len({p["table_sha256"] for p in ps}) != 1:
            why.append(f"V6: {key}: the parts used different tables")
        rows = h26.merge_rows(ps)[0]
        theirs = h27_parts.get(key)
        if theirs:
            w3, dist["V3"][key] = h28.gate(key, rows, h26.merge_rows(theirs)[0], mine=CONTROL, theirs=H27_ARM,
                                           truths=TRUTHS)
            why += [w.replace("V4: ", "V3: ", 1) for w in w3]
        else:
            why.append(f"V3: {key}: no committed row H27 part")
        mine28 = h28_parts.get(key)
        if mine28:
            why += u_scores(key, rows, h26.merge_rows(mine28)[0])
        else:
            why.append(f"V4: {key}: no committed row H28 part")
        why += identities(key, rows)
    return why, dist


# ------------------------------------------------------------------------------------------------ the tables
UNITS = {**{k: (k,) for k in TASK_KEYS}, "fold E10.5": TASK_KEYS}
PAIRS = ((f"{G} - {CONTROL}", G, CONTROL), (f"{U} - {CONTROL}", U, CONTROL), (f"{G} - {U}", G, U),
         (f"{G} - mean R", G, RS), *((f"{G} - {r}", G, r) for r in RS), (f"{G} - mean P", G, PS),
         *((f"{G} - {p}", G, p) for p in PS), (f"{G} - {F}", G, F), (f"{F} - {CONTROL}", F, CONTROL),
         (f"mean R - {CONTROL}", RS, CONTROL), (f"mean P - {CONTROL}", PS, CONTROL))


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: {lab: h30.compare(rows, ks, ts, a, b) for lab, a, b in PAIRS} for u, ks in UNITS.items()}
            for w, ts in WORLDS}


def statements(tab: dict) -> dict:
    """S1 to S6 on fold E10.5, per world, and whether each holds in both."""
    per = {}
    for w, _ in WORLDS:
        c = tab[w]["fold E10.5"]
        gc, gr, gp, gu, gf = (c[f"{G} - {CONTROL}"], c[f"{G} - mean R"], c[f"{G} - mean P"], c[f"{G} - {U}"],
                              c[f"{G} - {F}"])
        margin = gr["score"]["mean"]
        share = sum(gr[f"points:{g}"] for g in DE_GROUPS)
        per[w] = {
            "S1": bool(gc["score"]["mean"] >= BAR and gc["score"]["ahead"] >= MIN_AHEAD),
            "S2": bool(margin >= BAR and gr["score"]["ahead"] >= MIN_AHEAD
                       and all(c[f"{G} - {r}"]["score"]["mean"] > 0 for r in RS)),
            "S3": bool(gp["score"]["mean"] >= BAR and gp["score"]["ahead"] >= MIN_AHEAD),
            "S4": bool(all(gr[m] > 0 for m in DE_TERMS) and margin > 0 and share >= 0.5 * margin),
            "S5": bool(gu["score"]["mean"] >= 0),
            "S6": bool(gf["score"]["mean"] > 0 and gf["score"]["ahead"] >= F_MIN_AHEAD),
            "S4_margin": margin, "S4_de_groups_points": share}
    both = {s: all(per[w][s] for w, _ in WORLDS) for s in ("S1", "S2", "S3", "S4", "S5", "S6")}
    return {"per_world": per, "both": both, "all": bool(all(both.values()))}


def outcome(st: dict | None, void: bool) -> str:
    if void or st is None:
        return "VOID"
    return "AHEAD" if st["all"] else "NOT AHEAD"


FOLLOWS = {
    "VOID": "VOID: {why}. No figure exists under row H32 and there is no rerun under it; the heart file in force stands.",
    "AHEAD": ("AHEAD: S1 to S6 hold in both worlds on fold E10.5 (one truth; an in-sample diagnostic). What follows: "
              "one priced question to the captain with the file in force as the default: gain and loss per term, the "
              "practice value as a range from 0 (J6: test budget over practice 1.33 and 4), both readings if mixed, "
              "the Feng-filter doubt and the somite screen's result written out; nothing is built before the "
              "captain's A15 ruling and a yes by name."),
    "NOT AHEAD": ("NOT AHEAD: the fixed statements S1 to S6 do not all hold in both worlds. What follows: a one-line "
                  "note to the captain; the line is closed; no rerun under another row; the heart file in force "
                  "stands."),
}
L0 = ("Try H32 (run commit {commit}). A labelled in-sample diagnostic on one truth (E10.5, two backtests merged seed by "
      "seed); every figure is a practice reading on released stages and none is to be read beside a test score. "
      "Nothing is built and no file, rule, scope, size, count or depth changes. {bypass}.")
J5 = "J5: seed se about 0.01 to 0.04 on one truth (row H28), so the read is a labelled diagnostic."
J6 = ("J6: the trend's budget at the test ({t:g} day) over the practice budgets ({b}) is 1.33 and 4: any gain is shown "
      "as a range from zero to its practice value.")


def follows(res: dict) -> str:
    o = res["outcome"]
    if o == "VOID":
        return FOLLOWS["VOID"].format(why="; ".join(res["why"]))
    return FOLLOWS[o]


def reading(parts, h27_parts, h28_parts, sha: dict, commit: str | None, gate_ok: bool, g0: dict | None = None) -> dict:
    try:
        why, dist = void_reasons(parts, h27_parts, h28_parts, sha, commit, gate_ok)
    except (KeyError, TypeError, ValueError, IndexError) as e:  # a malformed part
        why, dist = [f"a part cannot be read ({type(e).__name__}: {e})"], {}
    n56 = n56_outcome()
    if why:
        res = {"void": True, "why": why, "guards": dist, "commit": commit, "outcome": "VOID",
               "n56_outcome_reported": n56}
        res["follows"] = follows(res)
        return res
    rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
    tab = tables(rows)
    st = statements(tab)
    res = {"void": False, "why": [], "row": ROW, "commit": commit, "guards": dist, "tables": tab, "statements": st,
           "outcome": outcome(st, False), "window_guard": WINDOW_BYPASS, "n56_outcome_reported": n56,
           "g0": g0, "cases": {k: [r["checks"]["cases"] for r in rows[k]] for k in TASK_KEYS},
           "sets_seed0": {k: rows[k][0]["checks"]["sets"] for k in TASK_KEYS}}
    res["follows"] = follows(res)
    return res


def _st(s: dict) -> str:
    return f"{s['mean']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10 (se {s['se']:.3f})"


def text(res: dict) -> str:
    if res["void"]:
        return res["follows"] + "\n" + f"Row N56's committed outcome (decides nothing): {res['n56_outcome_reported']}\n"
    lines = [L0.format(commit=res["commit"], bypass=WINDOW_BYPASS)]
    for key, g in res["guards"]["V3"].items():
        lines.append(f"GATE V3 {key} (the control against row H27's committed {H27_ARM}): {g['mean']:.3f} (ten-seed "
                     f"mean) / {g['seed']:.3f} (one seed); void beyond {GUARD_MEAN} / {GUARD_SEED}")
    for w, v in res["statements"]["per_world"].items():
        lines.append(f"STATEMENTS {w}: " + ", ".join(f"{s} {v[s]}" for s in ("S1", "S2", "S3", "S4", "S5", "S6"))
                     + f"; S4 margin {v['S4_margin']:+.3f}, in expression_change and cell_state "
                     f"{v['S4_de_groups_points']:+.3f}")
    lines += [f"OUTCOME: {res['outcome']}", res["follows"], J5,
              J6.format(t=TEST_BUDGET, b=", ".join(f"{k} {v:g} day" for k, v in PRACTICE_BUDGET.items())),
              f"Row N56's committed outcome, reported beside the reading (decides nothing): "
              f"{res['n56_outcome_reported']}."]
    if res.get("g0"):
        lines.append("G0 (committed): " + "; ".join(_g0_line("test pair", r) for r in res["g0"]["test_pair"].values()))
    lines.append("TABLES (decide nothing; paired by seed, mean [min, max] seeds ahead of 10 (se); per-term mean skill "
                 "differences; per-group board points):")
    for w, _ in WORLDS:
        for u in UNITS:
            for lab, c in res["tables"][w][u].items():
                lines.append(f"{TAG} {w} {u} {lab}: {_st(c['score'])}; "
                             + ", ".join(f"{m} {c[m]:+.4f}" for m in (*TERMS, *SHAPE_TERMS)) + "; "
                             + ", ".join(f"{g} {c[f'points:{g}']:+.3f}" for g in h30.GROUP_WEIGHT))
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def cmd_check(a) -> None:
    refuse_rows_and_register()
    print("rows and register: ok")


def cmd_tables(a) -> None:
    run_tables(Path(a.out), Path(a.data), Path(a.external), Path(a.h30) if a.h30 else None)


def cmd_gate(a) -> None:
    res = run_gate(Path(a.out), Path(a.tables), Path(a.data))
    sys.exit(0 if res["outcome"] == "PASS" else 3)


def cmd_run(a) -> None:
    refuse_run(Path(a.out) / PART)
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS}")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data), Path(a.tables))


def cmd_read(a) -> None:
    refuse_rows_and_register()
    out = Path(a.out)
    if (out / "reading.json").exists() or (out / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")
    if not COMMITTED_SHA.exists():
        sys.exit(f"{ROW} is refused: the tables' sha256 is not committed")
    gate_ok = committed_gate_ok(a.commit)
    g0 = json.loads(COMMITTED_GATE.read_text()) if COMMITTED_GATE.exists() else None
    sha = json.loads(COMMITTED_SHA.read_text())
    res = reading(load_parts(Path(a.runs)), h28.load_h27(Path(a.h27)), h28.load_parts(Path(a.h28)), sha, a.commit,
                  gate_ok, g0)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "reading.txt").write_text(text(res))
    print(f"OUTCOME: {res['outcome']}")
    print(res["follows"])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    t = sub.add_parser("tables")
    t.add_argument("--out", required=True)
    t.add_argument("--h30", default=None, help="row H30's tables folder, reused only as h30_reusable allows")
    t.add_argument("--data", default=str(harness.DATA_DIR))
    t.add_argument("--external", default=str(h30.EXTERNAL))
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--tables", required=True, help="the tables step's output folder")
    g.add_argument("--data", default=str(harness.DATA_DIR))
    r = sub.add_parser("run")
    r.add_argument("--task", required=True, help=f"one of {TASK_KEYS}")
    r.add_argument("--seeds", required=True, help="a range such as 0-1")
    r.add_argument("--out", required=True)
    r.add_argument("--tables", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--h27", default=str(h28.H27_PARTS))
    d.add_argument("--h28", default=str(H28_PARTS))
    d.add_argument("--out", required=True)
    d.add_argument("--commit", required=True, help="the run commit (full hash) the parts must come from")
    a = ap.parse_args(argv)
    {"check": cmd_check, "tables": cmd_tables, "gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
