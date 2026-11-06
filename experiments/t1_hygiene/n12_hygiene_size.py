#!/usr/bin/env python
"""Row N12 of experiments/harness/TRIES.md (T1; a diagnostic, not a try): R2's gate, measured before any build. Per
Feng direction group, how much of the top 400 up and the top 400 down genes by M1's factor fall in each of R2's
modules and in their union.

Target-free. It reads Feng's clean statistics through experiments/t1_m1/feng_dirs.py, as build_t1 does (X002,
manifest-checked), the committed extracts of X007 and X008 (sha256-checked against configs/t1_hygiene.yaml), the
pre-registered held-flat list, and the gene names of the released T1/E9.5_RNA.h5ad (names-only: through h5py, only the
var index element is read; no matrix, layer, obs or uns element). It reads no E10.5 file (R2: none is read before the
build), no target, no score, no board and no built file, and it counts, prints and stores no T1 cell. Its gene
universe is E9.5's names; the build maps the union by name onto E10.5's genes and records any union gene E10.5 lacks.

M1's factor, as build_t1 applies it (vec.t1_scope.external_shift, no clip): for a group g,
    log f = 0.4 x Feng's E9.5 -> E14.5 log ratio of mean CP10k (pseudocount 0.5),
0 where the gene has no Feng value and on the pre-registered held-flat list; 0.4 = feng_dirs.bracket_fraction(10.5,
12.5). Before anything is written the factors are checked against external_shift's own report on one synthetic cell per
group. The top lists are the 400 largest positive and the 400 most negative log f, ties by symbol (Python string
order).

R2's gate, from the count of the union U (configs/t1_hygiene.yaml: D, N2, and A1 and N1 while their lines are on) in
each top-400 list of each group:
- more than 133 in any list of any group: "captain" (nothing is built until the captain rules);
- fewer than 10 in every list of every group: "not built" (a near no-op: no rebuild, no s2 re-ruling);
- otherwise "build" (T1's s1 and s2 are rebuilt and format-checked, and T1's s2 goes back to the captain by 24 Nov).

Reported, deciding nothing, per group and list: the counts of D, A1, N1 and N2, each against a decile-matched null
(random sets of the module's mapped size matching its decile composition; deciles, K2b's bins, of the group's Feng mean
CP10k averaged over E9.5 and E14.5 among the genes with log f other than 0; one generator, default_rng(20261106),
draws every null in turn: groups by name, modules in the order D, A1, N1, N2, Table 5, 1,000 draws each, each draw its
bins from the lowest up; a module's mapped genes at log f = 0 can enter no list and are matched by genes that enter
none); the same at 100 and 1,000; abs(log f)-weighted shares; per group, the share of the sum of abs(log f) that U
holds and the Spearman correlation of log f with and without the hygiene over the genes with a Feng value; van den
Brink's Table 5 in place of D (U_table5); the ribosomal-protein, OXPHOS-pattern and mt- genes in each list and how many
U holds; module sizes, unmapped source symbols, and by name the module genes in each top-400 list and the module genes
on the pre-registered held-flat list. Each module's list (its genes mapped to the T1 names) is recorded with its sha256
(docs/PLAN.md 4.1, R2's provenance rule).

Outputs, in experiments/t1_hygiene/results/n12_hygiene_size/ (no target-derived content): n12_hygiene_size.json,
n12_hygiene_size.txt and hygiene_genes.csv (U, which the build reads through hygiene.load_union once its sha256 is
pinned). Each holds the inputs' sha256, the commit, the clean-tree flag (vec.provenance's rule, the outputs exempt,
from the start of the run to the write) and the UTC time; the text and the JSON also hold the sha256 of the outputs
written before them, and each output gets a vec.provenance build stamp (<file>.stamp.json) with its own sha256.

N12 runs once (R2: its result is not revisited after its numbers are seen; a change is a new row). Before anything is
read it refuses while configs/t1_hygiene.yaml pins union_sha256, or while any of its outputs or their stamps exists in
the output folder (a run that stopped early leaves files to remove by hand), unless --rerun-row names the new
experiments/harness/TRIES.md row that registers the rerun; that row is then recorded in every output.

Usage (DGX, repo root; only after the captain's sign-off on R2's gate, decision log):
  OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
      uv run --frozen --no-sync python -I -B experiments/t1_hygiene/n12_hygiene_size.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import warnings
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

from vec import harness, provenance, t1_scope
from vec.boards import REPO_ROOT

HERE = Path(__file__).resolve().parent
for _d in (HERE, REPO_ROOT / "experiments" / "t1_m1", REPO_ROOT / "experiments" / "t1_direction",
           REPO_ROOT / "experiments" / "external"):
    sys.path.insert(0, str(_d))
import dirlib as dl  # noqa: E402
import feng_dirs  # noqa: E402
import filterio as fio  # noqa: E402
import hygiene  # noqa: E402

OUT = HERE / "results" / "n12_hygiene_size"
TRIES = REPO_ROOT / "experiments" / "harness" / "TRIES.md"
ROW_ID = re.compile(r"[A-Z][0-9]+[a-z]?")  # a TRIES row's ID (R2, N12, K2b, ...)
T1_FILE = "T1/E9.5_RNA.h5ad"  # its gene names only; no E10.5 file is read before the build (R2)
NAMES = {"union": "hygiene_genes.csv", "text": "n12_hygiene_size.txt", "json": "n12_hygiene_size.json"}  # write order
SIDES = ("up", "down")
NULL_MODULES = (*hygiene.MODULES, hygiene.TABLE5)  # the null's draw order
UNIONS = ("U", "U_table5")  # R2's union; the same with Table 5 in place of D (reported only)
SETS = (*NULL_MODULES, *UNIONS)
LABEL = ("R2's gate (target-free: no E10.5 file, no target, no score, no board, no built file, no T1 cell count)")


# ---------------------------------------------------------------------------------------------------------------------
# The core: pure functions of given inputs


def m1_log_factor(direction, fraction: float, held: np.ndarray) -> np.ndarray:
    """M1's log factor for one group as build_t1 applies it (t1_scope.external_shift with no clip): `fraction` x the
    direction, 0 where the gene has no Feng value (NaN) and on the held genes."""
    logf = fraction * np.nan_to_num(np.asarray(direction, dtype=np.float64), nan=0.0)
    logf[np.asarray(held, dtype=bool)] = 0.0
    return logf


def rankings(logf: np.ndarray, genes: Sequence[str]) -> dict[str, np.ndarray]:
    """Gene positions with log f > 0, largest first ("up"), and with log f < 0, most negative first ("down"); ties by
    symbol, then by position should a name repeat. A side's top k is its first k entries; 0 is on neither side."""
    up = sorted(np.flatnonzero(logf > 0).tolist(), key=lambda i: (-logf[i], genes[i], i))
    down = sorted(np.flatnonzero(logf < 0).tolist(), key=lambda i: (logf[i], genes[i], i))
    return {"up": np.asarray(up, dtype=np.int64), "down": np.asarray(down, dtype=np.int64)}


def gate(counts: Mapping[str, Mapping[str, int]], upper: int, floor: int) -> str:
    """R2's gate from the union's count in each top list ({group: {"up": n, "down": n}}): "captain" when any list of
    any group holds more than `upper`; "not built" when every list of every group holds fewer than `floor`; otherwise
    "build"."""
    values = [int(v) for lists in counts.values() for v in lists.values()]
    if not values:
        raise ValueError("no top list to decide on")
    if max(values) > upper:
        return "captain"
    if max(values) < floor:
        return "not built"
    return "build"


def decile_bins(basis, universe: np.ndarray, bins: int) -> list[np.ndarray]:
    """K2b's bins (vec.harness.t3_decile_permutation): the universe's genes ranked by `basis` with a stable sort (ties
    keep gene order), cut into `bins` runs of consecutive ranks whose sizes differ by at most one, lowest first."""
    u = np.flatnonzero(universe)
    return np.array_split(u[np.argsort(np.asarray(basis, dtype=np.float64)[u], kind="stable")], bins)


def null_counts(members: np.ndarray, bins: Sequence[np.ndarray], lists: Mapping[str, np.ndarray], n_genes: int,
                draws: int, rng: np.random.Generator) -> tuple[list[int], dict[str, np.ndarray]]:
    """A module's decile-matched null: per draw, from each bin d in turn, k_d of its genes without replacement, k_d the
    module's genes in bin d. Returns k and, per list, the number of drawn genes in it at each draw."""
    bin_of = np.full(n_genes, -1, dtype=np.int64)
    for d, b in enumerate(bins):
        bin_of[b] = d
    own = bin_of[np.asarray(members, dtype=np.int64)]
    k = np.bincount(own[own >= 0], minlength=len(bins))
    masks = {}
    for name, idx in lists.items():
        masks[name] = np.zeros(n_genes, dtype=bool)
        masks[name][idx] = True
    out = {name: np.zeros(draws, dtype=np.int64) for name in lists}
    for t in range(draws):
        picks = [rng.choice(bins[d], int(k[d]), replace=False) for d in range(len(bins)) if k[d]]
        pick = np.concatenate(picks) if picks else np.zeros(0, dtype=np.int64)
        for name, m in masks.items():
            out[name][t] = int(m[pick].sum())
    return [int(x) for x in k], out


def null_summary(observed: int, draws: np.ndarray) -> dict:
    """The observed count against its null draws: their mean, sd, 5th and 95th percentiles, and the share of draws at
    or above the observed count ((1 + those draws) / (1 + all draws))."""
    d = np.asarray(draws, dtype=np.float64)
    return {"observed": int(observed), "null_mean": float(d.mean()),
            "null_sd": float(d.std(ddof=1)) if len(d) > 1 else 0.0,
            "null_q05": float(np.percentile(d, 5)), "null_q95": float(np.percentile(d, 95)),
            "p_at_or_above": float((1 + (d >= observed).sum()) / (1 + len(d)))}


def spearman(a: np.ndarray, b: np.ndarray) -> float | None:
    """Spearman's rho of a and b (None when it is undefined, as for a constant input)."""
    if len(a) < 2:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = float(spearmanr(a, b).statistic)
    return r if np.isfinite(r) else None


def n12(genes: Sequence[str], directions: Mapping[str, np.ndarray], basis: Mapping[str, np.ndarray],
        held_flat: Iterable[str], modules: Mapping[str, Sequence[str]], cfg: Mapping) -> dict:
    """N12 on given inputs; it reads and writes nothing.

    genes: the T1 gene names. directions: group -> Feng's E9.5 -> E14.5 log ratio of mean CP10k on `genes` (NaN where
    unmapped), as feng_dirs.feng_directions gives it. basis: group -> Feng's mean CP10k averaged over E9.5 and E14.5 on
    `genes`, the null's decile basis. held_flat: the pre-registered held-flat list. modules: source symbols of D, A1,
    N1, N2 and Table 5 (hygiene.build_modules with table5). cfg: configs/t1_hygiene.yaml.
    """
    genes = [str(g) for g in genes]
    n = len(genes)
    fraction = float(cfg["direction"]["fraction"])
    size, upper, floor = (int(cfg["gate"][k]) for k in ("list_size", "upper", "floor"))
    sizes = sorted({size, *(int(k) for k in cfg["n12"]["reported_sizes"])})
    nul = cfg["n12"]["decile_null"]
    held_set = {str(g) for g in held_flat}
    held = np.array([g in held_set for g in genes], dtype=bool)
    mapped = hygiene.map_modules(modules, genes)
    in_u = hygiene.union_modules(cfg)
    members = {"U": hygiene.union(mapped, in_u),
               "U_table5": hygiene.union(mapped, [hygiene.TABLE5 if m == "D" else m for m in in_u])}
    names = {m: set(mapped[m]["mapped"]) for m in NULL_MODULES} | {u: set(members[u]) for u in UNIONS}
    sets = {s: np.array([g in names[s] for g in genes], dtype=bool) for s in SETS}
    patterns = {p: np.array([re.search(rx, g) is not None for g in genes], dtype=bool)
                for p, rx in cfg["n12"]["reported_patterns"].items()}
    no_value = np.ones(n, dtype=bool)  # no Feng value in any group (unmapped in Feng)
    for d in directions.values():
        no_value &= np.isnan(np.asarray(d, dtype=np.float64))
    valueless = {genes[i] for i in np.flatnonzero(no_value)}
    rng = np.random.default_rng(int(nul["seed"]))
    groups = {}
    for g in sorted(directions):
        d = np.asarray(directions[g], dtype=np.float64)
        logf = m1_log_factor(d, fraction, held)
        ranks = rankings(logf, genes)
        lists = {f"{side}_{k}": ranks[side][:k] for side in SIDES for k in sizes}
        universe = logf != 0
        absf = np.abs(logf)
        rows = {}
        for name, idx in lists.items():
            w = absf[idx]
            rows[name] = {"length": int(len(idx)),
                          "counts": {s: int(sets[s][idx].sum()) for s in SETS},
                          "abs_log_f_share": {s: float(w[sets[s][idx]].sum() / w.sum()) if w.sum() > 0 else 0.0
                                              for s in SETS},
                          "patterns": {p: {"genes": int(pm[idx].sum()), "held_by_U": int((pm & sets["U"])[idx].sum())}
                                       for p, pm in patterns.items()}}
            if name.endswith(f"_{size}"):  # by name, in rank order, in the registered top lists
                rows[name]["module_genes"] = {m: [genes[i] for i in idx if sets[m][i]] for m in NULL_MODULES}
                rows[name]["pattern_genes"] = {p: [genes[i] for i in idx if pm[i]] for p, pm in patterns.items()}
        bins = decile_bins(basis[g], universe, int(nul["bins"]))
        null = {}
        for m in NULL_MODULES:
            pos = np.flatnonzero(sets[m])
            k, drawn = null_counts(pos, bins, lists, n, int(nul["draws"]), rng)
            null[m] = {"mapped": int(len(pos)), "in_universe": int(sum(k)), "per_decile": k,
                       "lists": {name: null_summary(int(sets[m][idx].sum()), drawn[name])
                                 for name, idx in lists.items()}}
        has_value = ~np.isnan(d)
        groups[g] = {"genes_with_feng_value": int(has_value.sum()), "nonzero_log_f": int(universe.sum()),
                     "positive": int(len(ranks["up"])), "negative": int(len(ranks["down"])),
                     "max_abs_log_f": float(absf.max()) if n else 0.0,
                     "gate_counts": {side: rows[f"{side}_{size}"]["counts"]["U"] for side in SIDES},
                     "held_share": {u: float(absf[sets[u]].sum() / absf.sum()) if absf.sum() > 0 else 0.0
                                    for u in UNIONS},
                     "spearman_with_without": {u: spearman(logf[has_value], np.where(sets[u], 0.0, logf)[has_value])
                                               for u in UNIONS},
                     "lists": rows, "null": null}
    counts = {g: r["gate_counts"] for g, r in groups.items()}
    verdict = gate(counts, upper, floor)
    flat = [(g, side, c) for g, r in counts.items() for side, c in r.items()]
    table5 = {g: {side: r["lists"][f"{side}_{size}"]["counts"]["U_table5"] for side in SIDES}
              for g, r in groups.items()}
    return {
        "verdict": verdict,
        "gate": {"rule": f"U holds more than {upper} of the top {size} up or down of any group: captain; fewer than "
                         f"{floor} in every list of every group: not built; otherwise build",
                 "counts": counts, "max": max(c for *_, c in flat),
                 "lists_above_upper": [f"{g} {s}" for g, s, c in flat if c > upper],
                 "lists_at_or_above_floor": [f"{g} {s}" for g, s, c in flat if c >= floor],
                 "with_table5_in_place_of_D": {"counts": table5, "reading": gate(table5, upper, floor),
                                               "decides": False}},
        "settings": {"fraction": fraction, "list_size": size, "upper": upper, "floor": floor, "sizes": sizes,
                     "null": {k: int(nul[k]) for k in ("draws", "seed", "bins")}, "union_modules": in_u,
                     "lines": {m: bool(cfg["lines"][m]) for m in hygiene.FROM_PATTERNS},
                     "patterns": dict(cfg["patterns"]), "reported_patterns": dict(cfg["n12"]["reported_patterns"]),
                     "t1_genes": n, "held_flat_genes_in_t1": int(held.sum())},
        "modules": {m: {"source": len(mapped[m]["source"]), "mapped": len(mapped[m]["mapped"]),
                        "genes": mapped[m]["mapped"], "sha256": hygiene.list_sha256(mapped[m]["mapped"]),
                        "source_sha256": hygiene.list_sha256(mapped[m]["source"]),
                        "unmapped": mapped[m]["unmapped"],
                        "mapped_on_held_flat": [x for x in mapped[m]["mapped"] if x in held_set],
                        "mapped_without_feng_value": [x for x in mapped[m]["mapped"] if x in valueless],
                        "in_union": m in in_u}
                    for m in NULL_MODULES},
        "unions": {u: {"genes": len(members[u]),
                       "per_module": {m: sum(m in ms for ms in members[u].values()) for m in NULL_MODULES
                                      if any(m in ms for ms in members[u].values())},
                       "on_held_flat": [x for x in members[u] if x in held_set]}
                   for u in UNIONS},
        "union_members": members["U"],
        "groups": groups,
    }


# ---------------------------------------------------------------------------------------------------------------------
# Inputs


def t1_gene_names(path: Path) -> list[str]:
    """The released file's gene names and nothing else: through h5py, only the var index element (anndata's dataframe
    encoding names it in var's `_index` attribute) is read, with anndata.io.read_elem; no matrix, layer, obs, obsm or
    uns element is read and no cell is counted. A file not in that encoding is refused rather than read whole."""
    with h5py.File(path, "r") as f:
        var = f.get("var")
        if not isinstance(var, h5py.Group) or "_index" not in var.attrs:
            sys.exit(f"{path}: var is not in anndata's dataframe encoding (no _index): its names cannot be read alone")
        return [str(g) for g in ad.io.read_elem(var[var.attrs["_index"]])]


def feng_inputs(genes: Sequence[str]) -> dict:
    """Feng's per-group direction on the T1 genes through feng_dirs.feng_directions, as build_t1 computes it (the clean
    reduction's statistics, manifest-checked), and from the same file the null's decile basis (each group's mean CP10k
    at E9.5 and E14.5, averaged) and Feng's gene symbols (A1 and N1). The directions are recomputed here from the same
    statistics and must be identical, or the run stops. Feng's cell counts are not kept."""
    path = feng_dirs.clean_stats_path()
    directions, rec = feng_dirs.feng_directions(list(genes), stats_path=path)
    if fio.verify_manifest_output(path).get("sha256") != rec["sha256"]:
        sys.exit("Feng's statistics changed between two reads")
    st = dl.load_stats(path)
    cfg = feng_dirs.load_config()
    idx = dl.sym_index(st["genes"], list(genes))
    eps = cfg["metrics"]["pseudocount_cp10k"]
    basis = {}
    for g in directions:
        a, b = (dl.feng_sel(st, stage, cfg["external_labels"][g][0]) for stage in feng_dirs.BRACKET)
        if not np.array_equal(dl.on_t1(dl.lfc_ratio(a, b, eps), idx), directions[g], equal_nan=True):
            sys.exit(f"{g}: Feng's direction recomputed here differs from feng_dirs.feng_directions'")
        basis[g] = dl.on_t1((a["cp"] + b["cp"]) / 2.0, idx)
    record = {k: rec[k] for k in ("file", "sha256", "bracket", "pseudocount_cp10k", "min_cells_ext", "mapped_genes",
                                  "of_genes")}
    record["groups"] = {g: {"feng_labels": list(rec["groups"][g]["feng_labels"])} for g in sorted(rec["groups"])}
    return {"directions": directions, "basis": basis, "symbols": [str(s) for s in st["genes"]], "record": record}


def tries_rows(path: Path = TRIES) -> set[str]:
    """The row IDs in experiments/harness/TRIES.md's tables: each table row's first word, where it is a row ID."""
    words = (ln[2:].split("|", 1)[0].split() for ln in path.read_text(encoding="utf-8").splitlines() if ln[:2] == "| ")
    return {w[0] for w in words if w and ROW_ID.fullmatch(w[0])}


def check_first_run(cfg: Mapping, outputs: Iterable[Path], rerun_row: str | None) -> None:
    """R2: N12's result is not revisited after its numbers are seen, and a change is a new row. Without a rerun row,
    refused while configs/t1_hygiene.yaml pins union_sha256 or while any of `outputs` exists; a rerun row must be a row
    of experiments/harness/TRIES.md other than N12 (registered before the rerun). Reads only the config's value, the
    outputs' existence and, for a rerun row, TRIES.md."""
    if rerun_row is not None:
        if rerun_row == "N12" or not ROW_ID.fullmatch(rerun_row) or rerun_row not in tries_rows():
            sys.exit(f"--rerun-row {rerun_row}: not a row of experiments/harness/TRIES.md other than N12 (a rerun is "
                     "registered there before it runs)")
        return
    why = "N12's result is not revisited (R2); a rerun is a new experiments/harness/TRIES.md row, given as --rerun-row"
    if cfg.get("union_sha256"):
        sys.exit(f"configs/t1_hygiene.yaml pins union_sha256: {why}")
    there = [p.name for p in outputs if Path(p).exists()]
    if there:
        sys.exit(f"N12's outputs exist already ({', '.join(there)}): {why} (the files of a run that stopped early are "
                 "removed by hand)")


def check_settings(cfg: Mapping) -> None:
    """The config's bracket fraction and pseudocount must be the ones feng_dirs uses."""
    d = cfg["direction"]
    want = feng_dirs.bracket_fraction(float(d["from_day"]), float(d["to_day"]))
    if not np.isclose(float(d["fraction"]), want, rtol=0, atol=1e-12):
        sys.exit(f"configs/t1_hygiene.yaml: fraction {d['fraction']} is not bracket_fraction = {want}")
    eps = feng_dirs.load_config()["metrics"]["pseudocount_cp10k"]
    if float(d["pseudocount_cp10k"]) != float(eps):
        sys.exit(f"configs/t1_hygiene.yaml: pseudocount {d['pseudocount_cp10k']} is not feng_dirs' {eps}")


def check_against_external_shift(genes: Sequence[str], directions: Mapping[str, np.ndarray], fraction: float,
                                 held_flat: Sequence[str], logfs: Mapping[str, np.ndarray]) -> None:
    """vec.t1_scope.external_shift itself, with build_t1's arguments, on one synthetic cell per group (a count of 1 on
    every gene; no data is read): its report's factor summary must match the log factors computed here, or the run
    stops (N12 would describe another factor than the build's)."""
    groups = sorted(directions)
    labels = [f"n12 {g}" for g in groups]
    C = sparse.csr_matrix(np.ones((len(groups), len(genes)), dtype=np.int32))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # repeated gene names, should any repeat, change nothing here
        cells = ad.AnnData(X=C.astype(np.float32), layers={"counts": C})
        cells.var_names = [str(g) for g in genes]
        cells.obs_names = [f"n12_{i}" for i in range(len(groups))]
        cells.obs[t1_scope.LABEL_KEY] = labels
        _, rep = t1_scope.external_shift(cells, directions, fraction, dict(zip(labels, groups, strict=True)),
                                         held_flat, seed=0)
    rows = {r["type"]: r for r in rep.rows if r["status"] == "shifted"}
    if set(rows) != set(logfs):
        sys.exit("external_shift moved other groups than the ones measured here")
    for g, r in rows.items():
        a = np.abs(logfs[g])
        same = ((r["n_up"], r["n_down"]) == (int((logfs[g] > np.log(1.1)).sum()), int((logfs[g] < -np.log(1.1)).sum()))
                and np.allclose([r["mean_abs_log_f"], r["max_abs_log_f"]], [a.mean(), a.max()], rtol=1e-12, atol=0))
        if not same:
            sys.exit(f"{g}: the log factor computed here differs from external_shift's report")


def tree_state(exempt: Sequence[Path]) -> dict:
    """HEAD and whether the tree is clean by vec.provenance's rule, the outputs exempt (None when git cannot tell)."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    return {"commit": commit, "clean": None if unclean is None else not unclean}


def _rel(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------------------------------------------------------
# The text report


def _null(s: dict) -> str:
    return (f"{s['observed']} (null {s['null_mean']:.1f}, 5-95% {s['null_q05']:.0f}-{s['null_q95']:.0f}, "
            f"p {s['p_at_or_above']:.3f})")


def _names(xs: Sequence[str]) -> str:
    return ", ".join(xs) if xs else "none"


def _r(x: float | None) -> str:
    return "n/a" if x is None else f"{x:.4f}"


def report_text(res: dict, inputs: dict, stamp: dict, outputs: Mapping[str, str]) -> str:
    st, gt = res["settings"], res["gate"]
    size = st["list_size"]
    why = {"captain": f"U holds more than {st['upper']} of a top-{size} list ({', '.join(gt['lists_above_upper'])}): "
                      "nothing is built until the captain rules",
           "not built": f"U holds fewer than {st['floor']} of {size} in every list of every group: the hygiene is a "
                        "near no-op and is not built (no rebuild, no s2 re-ruling)",
           "build": f"U holds {st['floor']} or more of {size} in some list and at most {st['upper']} in every list: "
                    "the hygiene is built, T1's s1 and s2 are rebuilt and format-checked, and T1's s2 goes back to "
                    "the captain by 24 Nov"}[res["verdict"]]
    lines = [f"N12 (experiments/harness/TRIES.md; T1, a diagnostic, not a try): {LABEL}",
             f"VERDICT: {res['verdict']} -- {why}",
             f"Gate: U = {' + '.join(st['union_modules'])} (A1 line {'on' if st['lines']['A1'] else 'off'}, N1 line "
             f"{'on' if st['lines']['N1'] else 'off'}), {res['unions']['U']['genes']:,} T1 genes; log f = "
             f"{st['fraction']} x Feng's E9.5 -> E14.5 log ratio of mean CP10k, 0 where unmapped and on the "
             f"{st['held_flat_genes_in_t1']:,} held-flat genes, no clip; top {size} up and down per group, ties by "
             "symbol"]
    lines += [f"  {g}: up {c['up']} of {res['groups'][g]['lists'][f'up_{size}']['length']}, down {c['down']} of "
              f"{res['groups'][g]['lists'][f'down_{size}']['length']}" for g, c in gt["counts"].items()]
    t5 = gt["with_table5_in_place_of_D"]
    lines.append(f"Reported only: with Table 5 in place of D the counts would read '{t5['reading']}' ("
                 + "; ".join(f"{g} up {c['up']}, down {c['down']}" for g, c in t5["counts"].items()) + ")")
    lines.append("Modules (source symbols -> mapped to the T1 names by exact symbol):")
    for m, info in res["modules"].items():
        what = "Table 5, in place of D (reported only)" if m == hygiene.TABLE5 else m
        lines.append(f"  {what}: {info['source']} -> {info['mapped']}{'' if info['in_union'] else ' (not in U)'}; "
                     f"unmapped: {_names(info['unmapped'])}; on the pre-registered held-flat list: "
                     f"{_names(info['mapped_on_held_flat'])}; mapped without a Feng value: "
                     f"{_names(info['mapped_without_feng_value'])}")
    lines.append("Per group (reported, deciding nothing; counts against the decile-matched null: observed (null mean, "
                 "5-95%, share of draws at or above)):")
    for g, r in res["groups"].items():
        lines.append(f"  {g}: {r['genes_with_feng_value']:,} genes with a Feng value, {r['nonzero_log_f']:,} with log "
                     f"f other than 0 ({r['positive']:,} up, {r['negative']:,} down), max |log f| "
                     f"{r['max_abs_log_f']:.3f}; U holds {r['held_share']['U']:.4f} of the sum of |log f| (Table 5 in "
                     f"place of D {r['held_share']['U_table5']:.4f}); Spearman of log f with and without the hygiene "
                     f"{_r(r['spearman_with_without']['U'])} ({_r(r['spearman_with_without']['U_table5'])})")
        for side in SIDES:
            for k in st["sizes"]:
                row, name = r["lists"][f"{side}_{k}"], f"{side}_{k}"
                nulls = "; ".join(f"{'Table 5' if m == hygiene.TABLE5 else m} {_null(r['null'][m]['lists'][name])}"
                                  for m in NULL_MODULES)
                shares = ", ".join(f"{s} {row['abs_log_f_share'][s]:.3f}" for s in ("U", *hygiene.MODULES))
                pats = ", ".join(f"{p} {v['genes']} ({v['held_by_U']} held)" for p, v in row["patterns"].items())
                lines.append(f"    top {k} {side} ({row['length']} genes): U {row['counts']['U']} (with Table 5 "
                             f"{row['counts']['U_table5']}); {nulls}; |log f| shares {shares}; {pats}")
            named = r["lists"][f"{side}_{size}"]
            lines.append(f"    module genes in the top {size} {side}: " + "; ".join(
                f"{'Table 5' if m == hygiene.TABLE5 else m}: {_names(xs)}" for m, xs in named["module_genes"].items()))
    lines.append("Outputs: " + "; ".join(f"{k} sha256 {v}" for k, v in outputs.items()))
    lines.append("Inputs: " + "; ".join(f"{k} {v['path']} sha256 {v['sha256']}" for k, v in inputs.items()))
    rerun = f", rerun under experiments/harness/TRIES.md row {stamp['rerun_row']}" if stamp.get("rerun_row") else ""
    lines.append(f"Commit {stamp['commit']}, tree clean {stamp['tree_clean']}, UTC {stamp['utc']}{rerun}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------------------------------------------------
# The run


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(OUT), help="output folder (default: the registered one)")
    ap.add_argument("--data", default=str(harness.DATA_DIR), help="the released data (default: data/raw)")
    ap.add_argument("--rerun-row", help="the experiments/harness/TRIES.md row that registers a rerun, recorded in "
                                        "every output; without it N12 refuses while configs/t1_hygiene.yaml pins "
                                        "union_sha256 or its outputs exist (R2: its result is not revisited)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    paths = {k: out / v for k, v in NAMES.items()}
    exempt = [*paths.values(), *(provenance.stamp_path(p) for p in paths.values())]
    cfg = hygiene.load_config()
    check_first_run(cfg, exempt, a.rerun_row)  # before anything else is read
    started = tree_state(exempt)
    check_settings(cfg)
    t1_path = Path(a.data) / T1_FILE
    genes = t1_gene_names(t1_path)
    feng = feng_inputs(genes)
    if not feng["directions"]:
        sys.exit("feng_dirs.feng_directions returned no group: there is no direction to measure")
    held_path = t1_scope.HELD_FLAT_GENES
    held_genes = t1_scope.held_flat_genes(held_path)
    fraction = float(cfg["direction"]["fraction"])
    held_set = set(held_genes)
    held = np.array([g in held_set for g in genes], dtype=bool)
    check_against_external_shift(genes, feng["directions"], fraction, held_genes,
                                 {g: m1_log_factor(d, fraction, held) for g, d in feng["directions"].items()})
    modules = hygiene.build_modules(feng["symbols"], cfg, table5=True)  # the extracts' sha256 checked
    res = n12(genes, feng["directions"], feng["basis"], held_genes, modules, cfg)
    inputs = {"config": {"path": _rel(hygiene.CONFIG), "sha256": hygiene.sha256_file(hygiene.CONFIG)},
              **{f"extract {m}": {"path": cfg["extracts"][m]["path"], "sha256": cfg["extracts"][m]["sha256"]}
                 for m in (*hygiene.FROM_EXTRACTS, hygiene.TABLE5)},
              "Feng statistics": {"path": feng["record"]["file"], "sha256": feng["record"]["sha256"]},
              "T1 gene names": {"path": T1_FILE, "sha256": hygiene.sha256_file(t1_path),
                                "names_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest()},
              "held-flat list": {"path": _rel(held_path), "sha256": hygiene.sha256_file(held_path)}}
    ended = tree_state(exempt)
    clean = (None if started["clean"] is None or ended["clean"] is None
             else started["clean"] and ended["clean"] and started["commit"] == ended["commit"])
    stamp = {"commit": ended["commit"], "tree_clean": clean, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             "rerun_row": a.rerun_row}
    rerun = f", rerun under experiments/harness/TRIES.md row {a.rerun_row}" if a.rerun_row else ""
    out.mkdir(parents=True, exist_ok=True)
    head = [f"N12 (experiments/harness/TRIES.md R2 and N12): R2's union U = "
            f"{' + '.join(res['settings']['union_modules'])}, {len(res['union_members'])} genes mapped to the T1 "
            "names, for the build's held list (experiments/t1_hygiene/hygiene.py load_union, once "
            "configs/t1_hygiene.yaml pins this file's sha256)",
            *(f"input {k}: {v['path']} sha256 {v['sha256']}" for k, v in inputs.items()),
            f"commit {stamp['commit']}, tree clean {stamp['tree_clean']}, utc {stamp['utc']}{rerun}"]
    written = {}
    data = hygiene.union_text(res["union_members"], head).encode()
    paths["union"].write_bytes(data)
    written[NAMES["union"]] = hashlib.sha256(data).hexdigest()
    text = report_text(res, inputs, stamp, written)
    paths["text"].write_bytes(text.encode())
    written[NAMES["text"]] = hashlib.sha256(text.encode()).hexdigest()
    body = {"label": f"N12 (experiments/harness/TRIES.md): {LABEL}", **res,
            "union_file": {"path": _rel(paths["union"]), "sha256": written[NAMES["union"]]},
            "feng": feng["record"], "inputs": inputs, "outputs": written, "provenance": stamp}
    paths["json"].write_bytes((json.dumps(body, indent=1) + "\n").encode())
    for p in paths.values():
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
