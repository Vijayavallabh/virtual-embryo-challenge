#!/usr/bin/env python
"""Try K26 (experiments/harness/TRIES.md; every step refuses without the row, on a dirty tree or at another HEAD): a
query-gene-conditioned response map. A rank-2 ridge map from a knocked-out gene's wild-type co-expression context (the
released T1 E8.5 and E9.5 files, experiments/t3_ext/context_map.py) to its frozen-form log response table, fitted over
the perturbations that have tables, evaluated at a fold's own knocked-out gene, multiplies the fold's base factor at
K24's target-free size exponent. A labelled reading for a ruling, not the adoption rule; A15 applies to any file.

Built on context_map.py, k18_built_files.py, k22_outside_forms.py (frozen-file check), k24_wider_tables.py (the tables'
loader, the size rule), outside_forms.py (product, permutation, figures), k21_factor_forms.py, k7_wt_revert.py and the
harness, all imported unchanged. No existing candidate, option or default changes.

Steps (each with --commit, the full hash of the run commit, which must be HEAD; outputs in the run folder R):
  features --out R    TARGET-FREE: the panel genes (WT E8.75), the label gene names (the pinned sc3 and bulk
                      manifests' units), T1 detection deciles, stage 0's random queries (20 per outside gene, seed
                      STAGE0_SEED), the q_rand draws (one per fold query, seed Q_RAND_SEED), every needed context;
                      R/features.json
  stage0 --out R      COUNTED DIAGNOSTIC (outside perturbations only; Mab21l2 and Gata4 excluded): leave-one-out, the
                      map against the LOO mean table and against its own random queries; R/stage0.json. Exit 7 when the
                      gate fails (stop: no fold read). Never recomputed over an existing stage0.json (its verdict is
                      read again)
  pregate --out R     TARGET-FREE (no fold's predicted table is compared with a table built from that fold's own target
                      file), after a passing stage 0: each fold's predicted table (Gata4 for g1/g2, fitted without
                      Gata4; Mab21l2 for x2, fitted without Mab21l2), sized to S* by K24's exponent, must change >= 100
                      genes and correlate < 0.9 with ln two, ln sc3, ln bulk and their mean (outside knockouts only)
                      and, on the g fold only, with ln of the barred factor. The barred factor (phi_M) is built from
                      T3/E9.5_mab21l2_ko.h5ad, the file x2 is scored against, so x2's pre-gate leaves it out
                      (PREGATE_TABLES). Else the arm and its controls drop (exit 5, nothing runs); R/pregate.json.
                      Each fold's map is trained with the OTHER fold's released knockout table as one label (the
                      specification's "fitted without" the fold's own gene), and a drop in either fold drops the arm
                      everywhere: that cross-fold dependence is stated in the row, not hidden
  run --part guard|score --world plain --task g1|g2 [--seeds] | --world dm [--workers --threads] | --world x2
                      the guard part (base arms only: barred_form in the g-cells, xko_base on x2) or a scoring part
                      (base, q_map, q_rand, q_perm); a scoring part refuses without a passing guard.json of this commit
  guard --runs R      the reproduction gate: barred_form against K18's committed rows in the four g-cells, xko_base
                      against K9's committed x2 rows of `raw` (try K2a's build), K17's tolerance; exit 6 when it fails
  read --runs R --out DIR   the fixed reading (never over an existing one: exit 4)

Exit codes: 0 done; 1 a crash (masked traceback); 2 refused; 3 INCOMPLETE (a file missing; nothing written); 4 a file
that is never overwritten exists; 5 pregate: the arm is dropped, nothing runs; 6 guard: the reproduction gate failed
(VOID); 7 stage 0's gate failed (stop). The script prints step lines only (no figure, no count): every figure goes to
the run folder's JSON files and the reading's files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
import traceback
from contextlib import contextmanager
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
T3X = REPO / "experiments" / "t3_ext"
for _p in (HERE, T3X):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
import context_map as cmap  # noqa: E402  (the features, the map, stage 0's statistics)
import k7_wt_revert as k7  # noqa: E402  (_paired, MIN_GAIN, MIN_AHEAD; the plain fold's preparation)
import k18_built_files as k18  # noqa: E402  (the worlds, tasks, stream, build path, scoring, guards; unchanged)
import k21_factor_forms as k21  # noqa: E402  (seed parser, K18's committed rows, _gated, _read_json; unchanged)
import k22_outside_forms as k22  # noqa: E402  (the frozen-file check, _stamp_ok; unchanged)
import k24_wider_tables as k24  # noqa: E402  (the tables' loader, folder checks, S* and the exponent; unchanged)
import outside_forms as ext  # noqa: E402  (product, permutation, figures, corr_between; unchanged)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

SCRIPT = "experiments/t3_decomp/k26_query_map.py"
ROW_SCRIPT = f"`{SCRIPT}`"
ROW_LABEL = "K26"
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"
ENV_PROCESSED = "VEC_T3_EXT_PROCESSED"  # the outside pipeline's processed/ folder (the per-gene unit files)
ENV_REFSTATS = "VEC_T3_EXT_REFSTATS"  # the reference statistics file the units were computed against

STAGE0_SEED = 20261122  # stage 0's random queries (the row)
Q_RAND_SEED = 20261126  # q_rand's draw, one gene per fold query (the row)
PERM_SEED = ext.PERM_SEED  # 20261113, q_perm's permutation within reference mean-expression deciles
QUERY = {"g": "Gata4", "x2": "Mab21l2"}  # each fold's knocked-out gene: g1/g2 (both worlds) and x2
FOLDS = ("g", "x2")
MIN_MOVED = 100  # the pre-gate: the sized table changes at least this many genes
MAX_CORR = 0.9  # the pre-gate: corr with each closed table, their mean and (g fold only) the barred factor below this
CLOSED = ("two", "sc3", "bulk")
# The tables each fold's pre-gate compares its predicted table with. two, sc3 and bulk are outside knockouts only (no
# fold's target enters them). The barred factor is built from T3/E9.5_mab21l2_ko.h5ad, x2's target file: x2 leaves it
# out, so the step stays target-free (row K26).
PREGATE_TABLES = {"g": (*CLOSED, "mean", "barred"), "x2": (*CLOSED, "mean")}
PREGATE_LEFT_OUT = {"g": {}, "x2": {"barred": "built from T3/E9.5_mab21l2_ko.h5ad, the file x2 is scored against"}}

STREAM = k18.STREAM
BASE = k18.BARRED  # "barred_form": the g-cells' base
BASE_X = "xko_base"  # x2's base: harness._t3_factor_xko, try K2a's build
ARM, RAND, PERM = "q_map", "q_rand", "q_perm"
TRY_ARMS = (ARM, RAND, PERM)  # the three counted tries
CONTROLS = (RAND, PERM)
G_ARMS = (BASE, *TRY_ARMS)
X_ARMS = (BASE_X, *TRY_ARMS)
TASKS = k18.TASKS  # g1, g2
WORLDS = k18.WORLDS  # plain, dm
CELLS = (("plain", "g1"), ("plain", "g2"), ("dm", "g1"), ("dm", "g2"), ("x2", "x2"))
X2_KEY = "T3/Gata4->Mab21l2"
DE_TERMS = ("de_score", "de_direction")
MIN_GAIN, MIN_AHEAD = k7.MIN_GAIN, k7.MIN_AHEAD  # 0.3 and 8 of 10
MAX_EXTRA_GATED = 1
MAX_TRIES = 3
VOID_TRIES = MAX_TRIES  # row K26: "a void counts": every VOID reading counts the three tries, whatever voided it
K18_SHA256 = dict(k21.K18_SHA256)
K9_X2 = HERE / "results" / "k9_residual" / "tasks" / "x2.json"
K9_X2_SHA256 = "8be12f798aeab27e68c44fd54fe3106ee949d8d2b4e4694bd9b5db0bfb5458ba"
K9_X2_ARM = "raw"
RESULTS = HERE / "results" / "k26_query_map"
EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_DROPPED, EXIT_GUARD, EXIT_STAGE0 = 2, 3, 4, 5, 6, 7
GUARD_FILE, GUARD_DIR = "guard.json", "guard"
_K26: dict = {}  # the sized tables per fold, set from pregate.json before any arm builds


class Refused(Exception):
    """A refusal of the whole step (exit status 2)."""


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _cm():
    import t3x_common as cm

    return cm


def _tb():
    import t3x_tables as tb

    return tb


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row and the preconditions


def row_present(path: Path | None = None) -> str:
    """The id of the one TRIES.md row that names this script; it must be K26 (exit 2 otherwise)."""
    path = TRIES_MD if path is None else path
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.startswith("| ") and ROW_SCRIPT in ln]
    ids = sorted({ln[2:].split(" ", 1)[0] for ln in rows})
    if ids != [ROW_LABEL]:
        refuse(f"{'no' if not ids else 'another or more than one'} row of {Path(path).name} names {ROW_SCRIPT} "
               f"({ids}); row {ROW_LABEL} is committed before any step")
    return ids[0]


DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing written; INCOMPLETE, not VOID")


def head_ok(st: dict, commit: str) -> None:
    if not commit or st["commit"] != commit:
        refuse(f"HEAD is {st['commit']}, not the run commit {commit}")


def env_folder(var: str, env: dict | None = None) -> Path:
    """A folder named by an environment variable, through outside_forms.table_dir's checks (k24.folder): set,
    absolute, no symlink, a directory, outside every checkout."""
    try:
        return k24.folder(var, env)
    except k24.Refused as e:
        raise Refused(str(e)) from None


def env_file(var: str, env: dict | None = None) -> Path:
    env = os.environ if env is None else env
    v = env.get(var)
    if not v:
        raise Refused(f"the environment variable {var} is not set: refused")
    p = Path(v)
    if not p.is_absolute() or p.is_symlink() or not p.is_file():
        raise Refused(f"{var} must name an absolute path of a plain file (no symlink): refused")
    return p


def preconditions(commit: str) -> tuple[str, dict]:
    """Before anything loads: the row, a clean tree, HEAD the run commit, the frozen T3 files at their pins, K24's
    listing at its pins and every outside folder or file variable valid."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    head_ok(st, commit)
    k22.frozen_ok()
    try:
        k24.read_listing()
        for var in (k24.ENV_TABLES, k24.ENV_TWO, ENV_PROCESSED):
            env_folder(var)
        env_file(ENV_REFSTATS)
    except (Refused, k24.Refused) as e:
        refuse(str(e))
    return rid, st


# ---------------------------------------------------------------------------------------------------------------------
# The labels


def unit_name(acc: str, gene: str) -> str:
    """process.write_unit's and bulk.write_bulk_unit's file stem."""
    return f"{acc}__{gene}".replace("+", "_and_").replace("/", "_")


def _manifest(d: Path, name: str, want_sha: str) -> dict:
    p = Path(d) / name / f"{name}_manifest.json"
    if p.is_symlink() or not p.is_file():
        raise Refused(f"{name}'s manifest is missing or a symlink: refused")
    raw = p.read_bytes()
    if sha256_bytes(raw) != want_sha:
        raise Refused(f"{name}'s manifest is not the one the pinned table's record names: refused")
    return json.loads(raw)


def unit_tables(manifest: dict, family: str, processed: Path, ref: dict, genes: list[str], datasets: dict,
                never: frozenset) -> dict[str, list[np.ndarray]]:
    """{knocked-out gene: [ln of each unit's finished pooled table]} from a pinned manifest's units. Each unit file is
    read once and checked against the manifest's sha256; its accession must have status 'cleared' in datasets.yaml
    (a held_phenocopy or excluded_design set, or one in NEVER_FETCH, refuses the step), its modality the family's, its
    genes the barred factor's and its reference statistics the manifest's. The table is finished as the pipeline
    finishes a unit (t3x_tables.finish_pooled, or finish over the whole-embryo pseudo-bulk for a bulk unit)."""
    tb = _tb()
    out: dict[str, list[np.ndarray]] = {}
    units = manifest.get("unit_sha256")
    if not isinstance(units, dict) or not units:
        raise Refused(f"the {family} manifest lists no units: refused")
    for key in sorted(units):
        acc, _, gene = key.partition(":")
        ds = datasets.get(acc)
        if acc in never or ds is None or ds.get("status") != "cleared":
            raise Refused(f"unit {key}: its dataset's status is {None if ds is None else ds.get('status')!r}, not "
                          "'cleared' (held and excluded sets are refused by status): refused")
        mod = "bulk" if family == "bulk" else "single_cell"
        if (ds.get("modality") or "single_cell") != mod:
            raise Refused(f"unit {key}: its dataset is not of the {family} family: refused")
        p = Path(processed) / acc / f"{unit_name(acc, gene)}.unit.json"
        if p.is_symlink() or not p.is_file():
            raise Refused(f"unit {key}: its unit file is missing or a symlink: refused")
        raw = p.read_bytes()
        if sha256_bytes(raw) != units[key]:
            raise Refused(f"unit {key}: the unit file is not at the manifest's sha256: refused")
        un = json.loads(raw)
        if un.get("genes") != list(genes) or un.get("refstats_sha256") != manifest.get("refstats_sha256") \
                or un.get("accession") != acc or un.get("gene") != gene:
            raise Refused(f"unit {key}: another gene list, reference or identity than the manifest's: refused")
        Lg = {g: np.array([np.nan if v is None else v for v in vals], dtype=np.float64)
              for g, vals in un["L"].items()}
        if family == "bulk":
            phi = tb.finish(Lg["__whole__"], tb.whole_den(ref))
        else:
            phi, _ = tb.finish_pooled(Lg, ref)
        out.setdefault(gene, []).append(np.log(np.asarray(phi, dtype=np.float64)))
    return out


def outside_labels(genes: list[str], env: dict | None = None) -> tuple[dict[str, np.ndarray], dict]:
    """The per-gene components behind sc3 and bulk as frozen-form log tables over `genes`: {gene: ln table}, the mean
    over a gene's units inside a family, then over the families that hold it. The tables are first loaded and checked
    through K24's loader at K24's pins (so the manifests are the pinned tables'); the reference statistics file must be
    the manifests' (sha256)."""
    cm = _cm()
    tb = _tb()
    try:
        loaded = k24.load_tables(list(genes), env=env)
    except k24.Refused as e:
        raise Refused(str(e)) from None
    processed = env_folder(ENV_PROCESSED, env)
    rs = env_file(ENV_REFSTATS, env)
    datasets = cm.load_datasets()
    never = frozenset(cm.NEVER_FETCH)
    per_family: dict[str, dict[str, list[np.ndarray]]] = {}
    rec: dict = {"manifests": {}, "units": {}}
    ref = None
    for key, family in (("sc3", "single_cell"), ("bulk", "bulk")):
        _, trec = loaded["tables"][key]
        d = env_folder(k24.SOURCES[key][0], env)
        man = _manifest(d, k24.SOURCES[key][1], trec["manifest_sha256"])
        if sha256_bytes(rs.read_bytes()) != man.get("refstats_sha256"):
            raise Refused(f"{ENV_REFSTATS} is not the reference statistics the {key} manifest names: refused")
        if ref is None:
            ref = tb.load_ref_stats(rs)
            if ref["genes"] != list(genes):
                raise Refused("the reference statistics' genes are not the barred factor's (order included): refused")
        per_family[family] = unit_tables(man, "bulk" if family == "bulk" else "sc", processed, ref, genes, datasets,
                                         never)
        rec["manifests"][key] = trec["manifest_sha256"]
        rec["units"][key] = sorted(man["unit_sha256"])
    out = {}
    for gene in sorted(set().union(*per_family.values())):
        fams = [np.mean(np.vstack(v[gene]), 0) for v in per_family.values() if gene in v]
        out[gene] = np.mean(np.vstack(fams), 0)
    rec["genes"] = sorted(out)
    rec["tables"] = {k: r["sha256"] for k, (_, r) in loaded["tables"].items()}
    rec["listing_sha256"] = loaded["listing_sha256"]
    return out, rec | {"loaded": loaded}


def frozen_labels(run_dir: Path, decision: Path, genes: list[str]) -> tuple[dict[str, np.ndarray], dict]:
    """Mab21l2's and Gata4's (pooled replicates) frozen-form tables from the pre-registered run's factor files through
    the builder's load_factor (sha256 against the decision file): {'Mab21l2': ln, 'Gata4': ln}."""
    bs = harness._t3_head()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    out, rec = {}, {}
    for gene, src in (("Mab21l2", "mab21l2"), ("Gata4", "gata4")):
        phi, r = bs.load_factor(Path(run_dir), bs.factor_name(src, dec["preset"]), list(genes), dec)
        out[gene] = np.log(np.asarray(phi, dtype=np.float64))
        rec[gene] = {k: r[k] for k in ("name", "sha256", "recorded_in")}
    return out, rec


def panel_genes(data: Path) -> list[str]:
    """The barred factor's gene order: WT E8.75 as the builder loads it (duplicates dropped)."""
    from vec.prep import load_stage

    return list(load_stage(Path(data) / "T2" / "E8.75.h5ad", drop_duplicates=True).var_names)


# ---------------------------------------------------------------------------------------------------------------------
# Features (target-free)


def plan_features(outside: list[str], deciles: dict[str, int]) -> dict:
    """Which genes get a context: the outside label genes with a context (in T1, at or above the floor, not a fold
    query), the two queries, stage 0's 20 random queries per outside gene (seed STAGE0_SEED, genes in sorted order, one
    generator) and one q_rand gene per fold query (seed Q_RAND_SEED, folds in FOLDS order). No draw is any label gene
    or query."""
    queries = [QUERY[f] for f in FOLDS]
    for q in queries:
        if q not in deciles:
            raise cmap.ContextRefused(f"the fold query {q} is absent from T1 or under the detection floor: refused")
    kept = sorted(g for g in outside if g in deciles and g not in queries)
    left = {g: ("a fold query" if g in queries else "absent from T1 or under the detection floor")
            for g in sorted(outside) if g not in kept}
    exclude = set(outside) | set(queries)
    rng = np.random.default_rng(STAGE0_SEED)
    s0 = {h: cmap.draw_same_decile(deciles, h, cmap.N_RANDOM, rng, exclude) for h in kept}
    rng2 = np.random.default_rng(Q_RAND_SEED)
    qr = {f: cmap.draw_same_decile(deciles, QUERY[f], 1, rng2, exclude)[0] for f in FOLDS}
    need = sorted(set(kept) | set(queries) | {g for v in s0.values() for g in v} | set(qr.values()))
    return {"outside_kept": kept, "outside_left_out": left, "stage0_random": s0, "q_rand": qr, "need": need}


def cmd_features(a) -> None:
    rid, st = preconditions(a.commit)
    out = Path(a.out)
    path = out / "features.json"
    if path.exists():  # never overwritten: a file of this commit stands (a rerun of the step), any other refuses
        if json.loads(path.read_text()).get("provenance", {}).get("commit") == a.commit:
            print("features: present", flush=True)
            sys.exit(0)
        print(f"refused: {path} exists from another commit (never overwritten)", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    out.mkdir(parents=True, exist_ok=True)
    try:
        genes = panel_genes(Path(a.data))
        try:
            labels, lrec = outside_labels(genes)
        except Refused as e:
            refuse(str(e))
        expr = cmap.load_t1(Path(a.data))
        det = cmap.detection(expr)
        level = cmap.mean_level(expr)
        deciles = cmap.decile_of(det, expr.genes)
        try:
            plan = plan_features(list(labels), deciles)
            ctx, crec = cmap.contexts(expr, plan["need"], genes, det=det, level=level)
        except cmap.ContextRefused as e:
            refuse(str(e))
        body = {"row": rid, "panel": genes, "plan": plan,
                "deciles": {g: deciles[g] for g in plan["need"]},
                "contexts": {g: [float(x) for x in v] for g, v in ctx.items()},
                "context_record": crec, "labels_record": {k: v for k, v in lrec.items() if k != "loaded"},
                "t1_files_sha256": {f: harness._t3_head().sha256(Path(a.data) / f) for f in cmap.T1_FILES},
                "constants": {"min_type_cells": cmap.MIN_TYPE_CELLS, "floor": cmap.DETECTION_FLOOR,
                              "rank": cmap.RANK, "lambda_per": cmap.LAMBDA_PER, "n_random": cmap.N_RANDOM,
                              "stage0_seed": STAGE0_SEED, "q_rand_seed": Q_RAND_SEED},
                "provenance": st}
        k18.write_task_file(path, body)
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print("features: written", flush=True)


def load_step_file(out: Path, name: str, commit: str) -> dict:
    p = Path(out) / name
    if not p.is_file():
        refuse(f"{p} is missing: the step before runs first")
    body = json.loads(p.read_text())
    if (body.get("provenance") or {}).get("commit") != commit or (body.get("provenance") or {}).get("tree_clean") \
            is not True:
        refuse(f"{name} is not of the run commit on a clean tree: refused")
    return body


# ---------------------------------------------------------------------------------------------------------------------
# Stage 0 (a counted diagnostic on outside perturbations only)


def stage0_from(feat: dict, labels: dict[str, np.ndarray]) -> dict:
    """Stage 0 from features.json's contexts and the outside labels (kept genes only)."""
    plan = feat["plan"]
    ctx = {g: np.asarray(v, dtype=np.float64) for g, v in feat["contexts"].items()}
    lab = {g: labels[g] for g in plan["outside_kept"]}
    return cmap.leave_one_out(ctx, lab, plan["stage0_random"], ctx)


def cmd_stage0(a) -> None:
    rid, st = preconditions(a.commit)
    out = Path(a.out)
    path = out / "stage0.json"
    if path.exists():  # never recomputed: the verdict of the existing file, if it is this commit's
        body = load_step_file(out, "stage0.json", a.commit)
        passes = (body.get("stage0") or {}).get("passes") is True
        print(f"stage 0: {'passes' if passes else 'fails'} (read again from {path.name})", flush=True)
        sys.exit(0 if passes else EXIT_STAGE0)
    feat = load_step_file(out, "features.json", a.commit)
    try:
        try:
            labels, lrec = outside_labels(feat["panel"])
        except Refused as e:
            refuse(str(e))
        if sorted(labels) != feat["labels_record"]["genes"]:
            refuse("the outside labels are not the ones features.json was made from: refused")
        res = stage0_from(feat, labels)
        body = {"row": rid, "stage0": res, "features_sha256": sha256_bytes((out / "features.json").read_bytes()),
                "labels_record": {k: v for k, v in lrec.items() if k != "loaded"}, "provenance": st}
        k18.write_task_file(path, body)
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print(f"stage 0: {'passes' if res['passes'] else 'fails'}", flush=True)
    if not res["passes"]:
        sys.exit(EXIT_STAGE0)


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate (target-free)


def sized(log_table) -> tuple[np.ndarray | None, float | None]:
    """exp(log_table) brought to K24's size S* by K24's exponent (k24.exponent); (None, None) when it has none."""
    P = np.exp(np.asarray(log_table, dtype=np.float64))
    e = k24.exponent(P)
    if e is None:
        return None, None
    return P ** e, float(e)


def fold_training(fold: str, outside: dict[str, np.ndarray], frozen: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """The fold's training perturbations: every outside gene and the released knockout that is NOT the fold's own."""
    other = {"g": "Mab21l2", "x2": "Gata4"}[fold]
    if QUERY[fold] in outside:
        raise Refused(f"an outside table holds the fold's own gene {QUERY[fold]}: refused")
    return {**outside, other: frozen[other]}


def pregate_fold(fold: str, ctx: dict[str, np.ndarray], train: dict[str, np.ndarray], rand_gene: str,
                 closed: dict[str, np.ndarray], phi_m: np.ndarray) -> dict:
    """One fold: the map fitted on `train`, evaluated at the fold's query (the arm) and at its q_rand gene, each sized;
    the drops (fewer than MIN_MOVED genes changed; corr >= MAX_CORR with a table of PREGATE_TABLES[fold]: the closed
    tables, their mean and, on the g fold only, phi_M; no exponent; the control unmakeable). phi_M is never read on
    x2's fold: it is built from x2's own target file."""
    if fold not in PREGATE_TABLES:
        raise RuntimeError(f"unknown fold {fold!r}")
    names = sorted(train)
    m = cmap.fit_map(np.vstack([ctx[g] for g in names]), np.vstack([train[g] for g in names]))
    arm, e = sized(m.predict(ctx[QUERY[fold]]))
    rnd, er = sized(m.predict(ctx[rand_gene]))
    why, corr = [], {}
    if arm is None:
        why.append("the predicted table has no size exponent (SD 0)")
    if rnd is None:
        why.append("q_rand's table has no size exponent (SD 0): the control cannot be made")
    figs = None
    if arm is not None:
        figs = ext.figures(np.ones_like(arm), arm)
        if figs["genes_ne_1"] < MIN_MOVED:
            why.append(f"the sized table changes fewer than {MIN_MOVED} genes")
        mean_closed = np.exp(np.mean(np.vstack([np.log(closed[k]) for k in CLOSED]), 0))
        avail = {**{k: closed[k] for k in CLOSED}, "mean": mean_closed}
        if "barred" in PREGATE_TABLES[fold]:
            avail["barred"] = phi_m
        for k in PREGATE_TABLES[fold]:
            c = ext.corr_between(arm, avail[k])
            corr[k] = c
            if c is None or c >= MAX_CORR:
                why.append(f"corr with {k} {'undefined' if c is None else 'at or above the bar'}")
    return {"train": names, "n_train": len(names), "rank": m.rank, "lambda": m.lam, "query": QUERY[fold],
            "rand_gene": rand_gene, "e": e, "e_rand": er, "figures": figs, "compared": list(PREGATE_TABLES[fold]),
            "left_out": dict(PREGATE_LEFT_OUT[fold]), "corr": corr, "dropped": why,
            "table": None if arm is None else [float(x) for x in arm],
            "table_rand": None if rnd is None else [float(x) for x in rnd]}


def pregate_all(ctx, outside, frozen, q_rand: dict[str, str], closed, phi_m) -> dict:
    folds = {f: pregate_fold(f, ctx, fold_training(f, outside, frozen), q_rand[f], closed, phi_m) for f in FOLDS}
    dropped = {f: v["dropped"] for f, v in folds.items() if v["dropped"]}
    # the arm is one try across all five cells: a drop in either fold drops it, with both controls, everywhere
    return {"folds": folds, "dropped": dropped, "runs": not dropped,
            "tries_if_run": len(TRY_ARMS) if not dropped else 0}


def cmd_pregate(a) -> None:
    rid, st = preconditions(a.commit)
    out = Path(a.out)
    path = out / "pregate.json"
    if path.exists():  # never overwritten: a file of this commit stands and its verdict is read again
        body = json.loads(path.read_text())
        if body.get("provenance", {}).get("commit") == a.commit:
            runs = (body.get("pregate") or {}).get("runs") is True
            print(f"pre-gate: present ({'the arm runs' if runs else 'the arm is dropped'})", flush=True)
            sys.exit(0 if runs else EXIT_DROPPED)
        print(f"refused: {path} exists from another commit (never overwritten)", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    feat = load_step_file(out, "features.json", a.commit)
    s0 = load_step_file(out, "stage0.json", a.commit)
    if (s0.get("stage0") or {}).get("passes") is not True:
        print("refused: stage 0 did not pass: no fold read", file=sys.stderr, flush=True)
        sys.exit(EXIT_STAGE0)
    if s0.get("features_sha256") != sha256_bytes((out / "features.json").read_bytes()):
        refuse("stage0.json was made from another features.json: refused")
    try:
        genes = feat["panel"]
        try:
            outside, lrec = outside_labels(genes)
        except Refused as e:
            refuse(str(e))
        outside = {g: outside[g] for g in feat["plan"]["outside_kept"]}
        frozen, frec = frozen_labels(Path(a.run_dir), Path(a.decision), genes)
        closed = {k: lrec["loaded"]["tables"][k][0] for k in CLOSED}
        k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
        if list(k18._BARRED["genes"]) != list(genes):
            raise RuntimeError("the barred factor's genes are not the panel's")
        ctx = {g: np.asarray(v, dtype=np.float64) for g, v in feat["contexts"].items()}
        try:
            pre = pregate_all(ctx, outside, frozen, feat["plan"]["q_rand"], closed, k18._BARRED["phi"])
        except Refused as e:
            refuse(str(e))
        body = {"row": rid, "pregate": pre, "panel": genes, "frozen_labels": frec,
                "features_sha256": s0["features_sha256"],
                "stage0_sha256": sha256_bytes((out / "stage0.json").read_bytes()),
                "tables": lrec["tables"], "listing_sha256": lrec["listing_sha256"], "provenance": st}
        k18.write_task_file(path, body)
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print(f"pre-gate: {'the arm runs' if pre['runs'] else 'the arm is dropped with its controls; nothing runs'}",
          flush=True)
    if not pre["runs"]:
        sys.exit(EXIT_DROPPED)


def load_pregate(out: Path, commit: str) -> dict:
    """pregate.json of this commit, the arm not dropped; fills _K26 with each fold's sized tables (panel order)."""
    body = load_step_file(out, "pregate.json", commit)
    pre = body["pregate"]
    if not pre.get("runs"):
        refuse("the pre-gate dropped the arm: nothing runs")
    _K26.clear()
    _K26.update(genes=list(body["panel"]), pregate_sha256=sha256_bytes((Path(out) / "pregate.json").read_bytes()),
                tables={f: {ARM: np.asarray(pre["folds"][f]["table"], dtype=np.float64),
                            RAND: np.asarray(pre["folds"][f]["table_rand"], dtype=np.float64)} for f in FOLDS})
    return body


# ---------------------------------------------------------------------------------------------------------------------
# The arms: the g-cells through K18's build path


def fold_tables(fold: str, ref_genes: list[str], basis) -> dict[str, np.ndarray]:
    """A fold's three try tables in the reference's gene order: q_map, q_rand, and q_perm (q_map permuted within
    deciles of `basis`, the world's reference mean expression, seed PERM_SEED)."""
    if not _K26.get("tables"):
        raise RuntimeError("the K26 tables are not loaded (pregate.json is read before any arm)")
    t = {a: k18.in_gene_order(_K26["tables"][fold][a], _K26["genes"], list(ref_genes)) for a in (ARM, RAND)}
    t[PERM] = ext.permuted(t[ARM], basis)
    return t


_K18_ORIGINAL = {"ARMS": k18.ARMS, "build_arm": k18.build_arm, "fold_factors": k18.fold_factors}


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K18's fold_factors (the original, unchanged) plus the g-fold's three tables in the reference's order."""
    F = _K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    F["k26"] = fold_tables("g", list(ref.genes), basis)
    return F


def arm_factor(arm: str, F: dict) -> np.ndarray:
    if arm == BASE:
        return F["phi_barred"]
    if arm not in TRY_ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    return ext.product(F["phi_barred"], 1.0, F["k26"][arm], 1.0)


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """As K24's build_arm: every arm draws from stable_rng(seed, STREAM) through bs.construct with the barred spec's
    keywords; only the factor passed differs; barred_form passes them unchanged."""
    if arm not in G_ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    if arm != BASE:
        kw["phi"] = arm_factor(arm, F)
    return bs.construct(ref, n, rng, **kw)


@contextmanager
def k18_arms(arms: tuple[str, ...]):
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's, restored on exit."""
    if not set(arms) <= set(G_ARMS) or BASE not in arms:
        raise RuntimeError(f"arms {arms} are not a set of this try's g-arms with the base")
    k18.ARMS, k18.build_arm, k18.fold_factors = tuple(arms), build_arm, fold_factors
    k18._FOLD.clear()
    try:
        yield
    finally:
        k18.ARMS = _K18_ORIGINAL["ARMS"]
        k18.build_arm = _K18_ORIGINAL["build_arm"]
        k18.fold_factors = _K18_ORIGINAL["fold_factors"]
        k18._FOLD.clear()


def g_candidates(arms: tuple[str, ...]) -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: k18._build_plain(ctx, arm),
                                 params={"form": arm, "stream": STREAM, "fold": "g", "query": QUERY["g"],
                                         "perm_seed": PERM_SEED if arm == PERM else None},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


# ---------------------------------------------------------------------------------------------------------------------
# The arms: x2 through try K2a's build


_X2: dict = {}


def x2_task():
    ts = [t for ts in harness.load_tasks().values() for t in ts if t.key == X2_KEY]
    if len(ts) != 1:
        raise RuntimeError(f"the harness holds {len(ts)} pseudo-tasks {X2_KEY!r}, not one")
    return ts[0]


def _x2_tables(ctx, prep: dict) -> dict:
    pool = ctx.inputs[ctx.task.reference_role]
    if _X2.get("prep") is not prep:
        _X2.clear()
        _X2.update(prep=prep, T=fold_tables("x2", list(prep["ref"].genes), harness.t3_mean_expression(pool)))
    return _X2["T"]


def build_x2(ctx, arm: str):
    """x2's arms. xko_base is harness._t3_factor_xko itself (try K2a's build). A try arm is the same build with K2a's
    factor times the arm's table: build_submission.construct on K2a's stream at the reference's own composition, then
    assemble and the harness's finish, as _t3_factor_xko calls them."""
    if arm == BASE_X:
        return harness._t3_factor_xko(ctx)
    if arm not in TRY_ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    bs = harness._t3_head()
    prep = harness._k2_prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    phi = ext.product(prep["phi"], 1.0, _x2_tables(ctx, prep)[arm], 1.0)
    built = bs.construct(ref, harness.count_for(ctx.board, pool.n_obs), bs.t3d.stable_rng(ctx.seed, harness.K2_STREAM),
                         f=ref.V.p, phi=phi)
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes),
                               bs.t3d.stable_rng(ctx.seed, f"{harness.K2_STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def x2_candidates(arms: tuple[str, ...]) -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: build_x2(ctx, arm),
                                 params={"form": arm, "stream": harness.K2_STREAM, "fold": "x2", "query": QUERY["x2"],
                                         "perm_seed": PERM_SEED if arm == PERM else None},
                                 requires=lambda t: t.key == X2_KEY)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


# ---------------------------------------------------------------------------------------------------------------------
# Running


def progress_log(part: str, world: str, task: str):
    lead = "K26 guard" if part == "guard" else "K26"

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"{lead} {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


def guard_passed(out: Path, st: dict) -> None:
    p = Path(out) / GUARD_FILE
    if not p.is_file():
        refuse(f"{p} is missing: the guard part runs and passes first")
    g = json.loads(p.read_text())
    prov = g.get("provenance") or {}
    if g.get("passes") is not True:
        refuse("the reproduction gate did not pass: the run is VOID and no scoring part runs")
    if prov.get("commit") != st["commit"] or prov.get("machine") != st["machine"]:
        refuse("guard.json is from another commit or machine: refused")


def _harness_part(a, st: dict, out: Path, part: str, world: str, task, cands, folder_: Path) -> None:
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    lead = "K26 guard" if part == "guard" else "K26"
    for s in k21.parse_seeds(a.seeds):
        path = folder_ / f"seed_{s}.json"
        if path.exists():
            print(f"{lead} {world} {a.task or world} seed {s}: present, skipped", flush=True)
            continue
        res = harness.run_task(task, cands, seeds=[s], data_dir=Path(a.data),
                               log=progress_log(part, world, a.task or world))
        got = res["provenance"]
        if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
            raise k22._Unclean(f"the harness's stamp says commit {got.get('commit')}, clean {got.get('tree_clean')}: "
                               "no task file written (INCOMPLETE)")
        k22._stamp_ok(st)
        res["provenance"]["machine"] = st["machine"]
        res |= {"world": world, "part": part, "task_source": SCRIPT, "seed": s,
                "pregate_sha256": _K26["pregate_sha256"]}
        k18.write_task_file(path, res)
        print(f"wrote {path}", flush=True)


def run_plain(a, st: dict, out: Path, part: str) -> None:
    root = out / GUARD_DIR if part == "guard" else out
    k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    if list(k18._BARRED["genes"]) != _K26["genes"]:
        raise RuntimeError("the barred factor's genes are not pregate.json's panel")
    arms = (BASE,) if part == "guard" else G_ARMS
    with k18_arms(arms):
        _harness_part(a, st, out, part, "plain", TASKS[a.task], g_candidates(arms), root / f"plain_{a.task}")


def run_x2(a, st: dict, out: Path, part: str) -> None:
    root = out / GUARD_DIR if part == "guard" else out
    arms = (BASE_X,) if part == "guard" else X_ARMS
    _X2.clear()
    _harness_part(a, st, out, part, "x2", x2_task(), x2_candidates(arms), root / "x2")


def run_dm(a, st: dict, out: Path, part: str) -> None:
    """The depth-matched world, set up as k24.run_dm (k23's set-up), each (replicate, seed) row through k18.dm_row
    with this part's arms; rows already written at this commit on this machine are skipped."""
    import multiprocessing as mp

    import orth_diag as od

    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / "dm"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    seeds = k21.parse_seeds(a.seeds)
    bs = harness._t3_head()
    t3d = bs.t3d
    spec, dec, rules = k18.barred_spec(Path(a.decision))
    if list(rules["protocol"]["seeds"]) != list(SEEDS):
        raise RuntimeError("the frozen script's seeds are not vec.adoption.SEEDS")
    paths = od.context_paths(Path(a.data), with_targets=True)
    n_par = a.workers * a.threads
    by_name = {v: k for k, v in k18.DM_TARGETS.items()}
    arms = (BASE,) if part == "guard" else G_ARMS
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=Path(tmp))
        if sorted(T.name for T in ctx.targets) != sorted(k18.DM_TARGETS.values()):
            raise RuntimeError("the depth-matched targets are not rep1 and rep2")
        k18.quiet_selfcheck(t3d, ctx)
        genes = list(ctx.genes)
        phi_b, frec = k18.barred_factor(Path(a.run_dir), dec, spec, genes)
        ref_C, ref_lab, ref_P, gap = od.reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        wt_path = Path(paths["wt"])
        series = bs.wt_revert_series(wt_path, {n: Path(a.data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                               if n != wt_path.stem})
        phi, per, group, rrec = bs.wt_revert_factors(wt_path, series, genes, ctx.wt.pool.C)
        refw = bs.Reference(preset=ctx.preset_name, rules=rules, genes=genes, V=ctx.V, pool=ctx.wt.pool,
                            coords=np.zeros((ctx.wt.pool.C.shape[0], 3)))
        files = [*(Path(a.data) / f for f in od.WT_FILES.values()), *(Path(a.data) / f for f in od.KO_FILES.values()),
                 *series.values(), Path(a.decision), Path(a.run_dir) / "factors" / f"{od.BARRED}.csv"]
        setup = {"world": "dm", "part": part, "barred_factor": frec,
                 "wt_revert": {k: v for k, v in rrec.items() if k != "series"},
                 "reference_counts_largest_gap": gap, "pregate_sha256": _K26["pregate_sha256"],
                 "files_sha256": {str(f): t3d.sha256(f) for f in dict.fromkeys(files)}}
        sp = folder_ / "setup.json"
        if sp.exists():
            old = json.loads(sp.read_text())
            if {k: old.get(k) for k in setup} != json.loads(json.dumps(setup, default=float)):
                sys.exit(f"{sp} records another set-up: not run")
        else:
            k18.write_task_file(sp, setup | {"provenance": st})
        with k18_arms(arms):
            F = k18.fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b,
                                 spec)
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder_ / f"{by_name[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                k22._stamp_ok(st)
                path = folder_ / f"{by_name[row['target']]}_seed_{row['seed']}.json"
                k18.write_task_file(path, row | {"world": "dm", "part": part, "task": by_name[row["target"]],
                                                 "arms": list(arms), "pregate_sha256": _K26["pregate_sha256"],
                                                 "provenance": st, "task_source": SCRIPT})
                print(f"wrote {path}", flush=True)

            if a.workers > 1 and len(todo) > 1:
                if not bs.openmp_limited():
                    raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
                with mp.get_context("fork").Pool(min(a.workers, len(todo))) as p:
                    for row in p.imap_unordered(k18._dm_task, todo):
                        write(row)
            else:
                for t in todo:
                    write(k18._dm_task(t))


def cmd_run(a) -> None:
    if a.world == "plain" and a.task not in TASKS:
        refuse("--world plain needs --task g1 or g2")
    if a.world != "plain" and a.task is not None:
        refuse("--task is for --world plain only")
    _, st = preconditions(a.commit)
    out = Path(a.out)
    s0 = load_step_file(out, "stage0.json", a.commit)
    if (s0.get("stage0") or {}).get("passes") is not True:
        refuse("stage 0 did not pass: no fold read")
    if a.part == "guard":
        try:
            guard_refs()
        except Refused as e:
            refuse(str(e))
    else:
        guard_passed(out, st)
    load_pregate(out, a.commit)
    try:
        {"plain": run_plain, "dm": run_dm, "x2": run_x2}[a.world](a, st, out, a.part)
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seeds written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)


# ---------------------------------------------------------------------------------------------------------------------
# The reproduction gate


def k9_rows(path: Path = K9_X2, pin: str = K9_X2_SHA256) -> list[dict]:
    """K9's committed x2 rows (try K2a's build under the name `raw`), at their pinned sha256."""
    p = Path(path)
    if not p.is_file() or sha256_bytes(p.read_bytes()) != pin:
        raise Refused(f"K9's committed {p.name} is missing or not at its pinned sha256: refused")
    return json.loads(p.read_text())["rows"]


def guard_refs(k18_dir: Path = k21.K18_TASKS, k9_path: Path = K9_X2) -> dict:
    saved = dict(k21.K18_SHA256)
    k21.K18_SHA256.clear()
    k21.K18_SHA256.update(K18_SHA256)
    try:
        k18r, why = k21.k18_rows(k18_dir)
    finally:
        k21.K18_SHA256.clear()
        k21.K18_SHA256.update(saved)
    if why:
        raise Refused("; ".join(why) + ": refused")
    return {"k18": k18r, "k9_x2": k9_rows(k9_path)}


def guard_files() -> list[str]:
    return ([f"{GUARD_DIR}/plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"{GUARD_DIR}/dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"{GUARD_DIR}/x2/seed_{s}.json" for s in SEEDS])


def score_files() -> list[str]:
    return ([f"plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"x2/seed_{s}.json" for s in SEEDS])


def base_of(cell: tuple[str, str]) -> str:
    return BASE_X if cell[0] == "x2" else BASE


def cell_of_file(rel: str) -> tuple[tuple[str, str], int]:
    """((world, task), seed) of a task file's path relative to its part's root."""
    seed = int(rel.rsplit("_", 1)[1][:-len(".json")])
    if rel.startswith("plain_"):
        return ("plain", rel.split("/")[0][len("plain_"):]), seed
    if rel.startswith("dm/"):
        return ("dm", rel[3:5]), seed
    if rel.startswith("x2/"):
        return ("x2", "x2"), seed
    raise ValueError(f"unknown task file {rel}")


def rows_of(runs: Path, files: list[str], prefix: str, want: dict, why: list[str], sha: dict, prov: dict,
            pregate_sha: str | None = None) -> dict:
    """Rows per cell from task files (plain and x2: rows[0]; dm: the row), each checked for its task, seed, arms and
    the pre-gate it was built under."""
    rows: dict = {c: [] for c in CELLS}
    x2_key = None
    for f in files:
        r, sha[f] = k21._read_json(Path(runs) / f)
        prov[f] = r["provenance"]
        cell, seed = cell_of_file(f[len(prefix):])
        if pregate_sha is not None and r.get("pregate_sha256") != pregate_sha:
            why.append(f"{f}: built under another pregate.json")
        if cell[0] == "dm":
            if r.get("task") != cell[1] or r.get("seed") != seed:
                why.append(f"{f}: another task or seed")
                continue
            row, need = r, want[cell]
        else:
            key = TASKS[cell[1]].key if cell[0] == "plain" else X2_KEY
            x2_key = key if cell[0] == "x2" else x2_key
            if r["task"]["key"] != key or [q["seed"] for q in r["rows"]] != [seed]:
                why.append(f"{f}: another task or seed")
                continue
            row, need = r["rows"][0], (harness.FLOOR, *want[cell])
        miss = [x for x in need if x not in row.get("scores", {}) or x not in row.get("skills", {})
                or x not in row.get("t3", {})]
        if miss:
            why.append(f"{f}: missing arms {miss}")
            continue
        rows[cell].append(row)
    return {k: sorted(v, key=lambda q: q["seed"]) for k, v in rows.items()}


def guard_eval(rows: dict, refs: dict) -> dict:
    """The reproduction gate: barred_form within K17's tolerance (k18.guard_passes: ten-seed mean within 0.1, every seed
    within 0.5) of K18's committed rows in the four g-cells, and xko_base within it of K9's committed `raw` on x2."""
    out = {}
    for c in CELLS:
        rr = rows.get(c) or []
        ok_seeds = [q["seed"] for q in rr] == list(SEEDS)
        if c[0] == "x2":
            d = k18.distance(rr, BASE_X, refs["k9_x2"], K9_X2_ARM) if ok_seeds else None
        else:
            d = k18.distance(rr, BASE, refs["k18"].get(c) or [], BASE) if ok_seeds else None
        out[f"{c[0]} {c[1]}"] = {"distance": d, "passes": bool(ok_seeds and k18.guard_passes(d))}
    failed = sorted({BASE_X if k.startswith("x2") else BASE for k, v in out.items() if not v["passes"]})
    return {"cells": out, "failed": failed, "passes": not failed}


def cmd_guard(a) -> None:
    rid, st = preconditions(a.commit)
    runs = Path(a.runs)
    p = runs / GUARD_FILE
    if p.exists():
        print(f"refused: {p} exists: the gate is made once", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    missing = [f for f in guard_files() if not (runs / f).exists()]
    if missing:
        print("INCOMPLETE: guard files missing", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    try:
        refs = guard_refs()
    except Refused as e:
        refuse(str(e))
    pre = load_step_file(runs, "pregate.json", a.commit)
    psha = sha256_bytes((runs / "pregate.json").read_bytes())
    try:
        why: list[str] = []
        sha, prov = {}, {}
        want = {c: (base_of(c),) for c in CELLS}
        rows = rows_of(runs, guard_files(), f"{GUARD_DIR}/", want, why, sha, prov, psha)
        bad = [f for f, q in prov.items() if q.get("tree_clean") is not True or q.get("commit") != st["commit"]
               or q.get("machine") != st["machine"]]
        if bad or why:
            refuse(f"the guard files are not this commit's and machine's, or not well formed ({(bad + why)[:4]})")
        res = guard_eval(rows, refs)
        k18.write_task_file(p, {**res, "row": rid, "files_sha256": sha, "pregate_sha256": psha,
                                "pregate_runs": pre["pregate"]["runs"], "provenance": st})
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print(f"reproduction gate: {'passes' if res['passes'] else 'FAILS'}", flush=True)
    if not res["passes"]:
        sys.exit(EXIT_GUARD)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def summarise(rows: list[dict], arms: tuple[str, ...], base: str) -> dict:
    """Per arm: mean and sd of the score, mean skill per metric, severity states, gated seeds; paired (k7._paired):
    each try arm against the base, q_map against each control, the base against the floor where present; per-term
    skill differences against the base; the de_score plus de_direction skill change of each try arm."""
    out: dict = {"arms": {}, "paired": {}, "per_term": {}, "de_change": {}}
    for x in arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "severity_states": {v: sum(z["state"] == v for z in sev) for v in sorted({z["state"] for z in sev})},
            "severity_gated": k21._gated(sev)}
    pairs = [(x, base) for x in TRY_ARMS if x in arms] + [(ARM, c) for c in CONTROLS if c in arms and ARM in arms]
    if harness.FLOOR in arms:
        pairs.append((base, harness.FLOOR))
    out["paired"] = {f"{x} - {y}": k7._paired(rows, x, y) for x, y in pairs}
    for x in TRY_ARMS:
        if x not in arms:
            continue
        out["per_term"][x] = {m: float(np.mean([q["skills"][x][m] - q["skills"][base][m] for q in rows]))
                              for m in rows[0]["skills"][base]}
        out["de_change"][x] = float(np.mean([sum(q["skills"][x][m] - q["skills"][base][m] for m in DE_TERMS)
                                             for q in rows]))
    return out


def _beats(v: dict | None) -> bool:
    return bool(v is not None and v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD)


def decide(per: dict) -> dict:
    """The pass line fixed in the row, over the five cells (per[(world, task)] = summarise()):
    (a) AHEAD (paired mean > 0.3 and >= 8 of 10 seeds ahead of the base) on x2 and on at least one g-cell in each
        world; (b) BEHIND (paired mean < -0.3) in none; (c) gated seeds of q_map at most the base's plus one in every
        cell; (d) q_map beats q_rand AND q_perm by > 0.3 with 8 of 10 in every AHEAD cell; (e) the de_score plus
        de_direction skill change of q_map against the base positive in every AHEAD cell."""
    key = lambda c: f"{c[0]} {c[1]}"  # noqa: E731
    d = {c: per[c]["paired"][f"{ARM} - {base_of(c)}"] for c in CELLS}
    ahead = {c: _beats(v) for c, v in d.items()}
    behind = {c: bool(v["mean"] < -MIN_GAIN) for c, v in d.items()}
    gated = {c: (per[c]["arms"][ARM]["severity_gated"], per[c]["arms"][base_of(c)]["severity_gated"]) for c in CELLS}
    up = [c for c in CELLS if ahead[c]]
    a = (ahead[("x2", "x2")] and any(ahead[("plain", n)] for n in TASKS) and any(ahead[("dm", n)] for n in TASKS))
    b = not any(behind.values())
    cc = all(g <= g0 + MAX_EXTRA_GATED for g, g0 in gated.values())
    dd = all(_beats(per[c]["paired"].get(f"{ARM} - {x}")) for c in up for x in CONTROLS)
    ee = all(per[c]["de_change"][ARM] > 0 for c in up)
    passes = bool(a and b and cc and dd and ee)
    return {"passes": passes, "a": bool(a), "b": bool(b), "c": bool(cc), "d": bool(dd), "e": bool(ee),
            "ahead_cells": [key(c) for c in up], "behind_cells": [key(c) for c in CELLS if behind[c]],
            "means": {key(c): float(v["mean"]) for c, v in d.items()},
            "gated_arm_base": {key(c): list(v) for c, v in gated.items()},
            "outcome": ("PASS: the kill test is passed by the line fixed in row K26; a labelled reading for a ruling, "
                        "not the adoption rule; A15 applies (a priced question to the captain); no file is built here"
                        if passes else "FAIL: the kill test is not passed; no file; T3's built and listed files are "
                        "unchanged")}


INCOMPLETE = ("not a reading: a step file or a seed's task file is missing (a crash, a refusal such as a dirty tree, "
              "out of memory or a kill); nothing was read. Rerun the missing part into this folder at the run's commit "
              "on the same machine (the runs skip what is there); after a committed code fix, rerun every step at the "
              "new commit into a fresh folder and never read this one")


def reading_code() -> dict:
    st = k18.stamp()
    h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
    return {"script": SCRIPT, "script_sha256": h(Path(__file__)), "context_map_sha256": h(T3X / "context_map.py"),
            "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}


def reading(runs: Path, refs: dict | None = None) -> dict:
    """Apply the reading of row K26 to a run folder: a STOP at stage 0 or at the pre-gate is a reading of its own (no
    fold read, no try beyond stage 0); otherwise the guard part is read again, then the scoring part."""
    runs = Path(runs)
    why: list[str] = []
    sha, prov, body = {}, {}, {}

    def load(*names) -> list[str]:
        miss = [f for f in names if not (runs / f).exists()]
        for f in names:
            if f not in miss and f not in body:
                body[f], sha[f] = k21._read_json(runs / f)
                prov[f] = body[f]["provenance"]
        return miss

    miss = load("features.json", "stage0.json")
    if miss:
        return {"incomplete": True, "void": False, "missing": miss, "outcome": INCOMPLETE}
    feat, s0 = body["features.json"], body["stage0.json"]
    s0res = s0.get("stage0") or {}
    fplan = feat.get("plan") or {}
    shown = (*QUERY.values(), *(fplan.get("q_rand") or {}).values())
    common = {"stage0": s0res, "breadth": {gn: v for gn, v in ((feat.get("context_record") or {}).get("breadth")
                                                               or {}).items() if gn in shown},
              "types_used": (feat.get("context_record") or {}).get("types"),
              "outside_left_out": fplan.get("outside_left_out")}
    if s0.get("features_sha256") != sha["features.json"]:
        why.append("stage0.json was made from another features.json")
    later = [f for f in ("pregate.json", GUARD_FILE) if (runs / f).exists()]
    if s0res.get("passes") is not True:
        stop = {"incomplete": False, **common, "files_sha256": dict(sha), "provenance": dict(prov),
                "tries_counted": 0, "stage0_counted": 1, "stopped": "stage0"}
        if later or why:
            return stop | {"void": True, "tries_counted": VOID_TRIES,
                           "why": [*why, *(f"{f} exists after a stop at stage 0" for f in later)],
                           "outcome": "void: the run went past the stop at stage 0, or its files disagree; the three "
                                      "tries count"}
        return stop | {"void": False, "why": [], "outcome": "STOP at stage 0: the gate failed; no fold read, no file"}
    miss = load("pregate.json")
    if miss:
        return {"incomplete": True, "void": False, "missing": miss, "outcome": INCOMPLETE}
    pre = body["pregate.json"]
    common["pregate"] = _pregate_report(pre.get("pregate") or {})
    if pre.get("features_sha256") != sha["features.json"]:
        why.append("pregate.json was made from another features.json")
    if pre.get("stage0_sha256") != sha["stage0.json"]:
        why.append("pregate.json was made from another stage0.json")
    if not (pre.get("pregate") or {}).get("runs"):
        stop = {"incomplete": False, **common, "files_sha256": dict(sha), "provenance": dict(prov),
                "tries_counted": 0, "stage0_counted": 1, "stopped": "pregate"}
        if (runs / GUARD_FILE).exists() or why:
            return stop | {"void": True, "tries_counted": VOID_TRIES,
                           "why": [*why] + (["guard.json exists after the pre-gate's drop"]
                                            if (runs / GUARD_FILE).exists() else []),
                           "outcome": "void: the run went past the pre-gate's drop, or its files disagree; the three "
                                      "tries count"}
        return stop | {"void": False, "why": [], "outcome": "STOP at the pre-gate: the arm is dropped with its "
                                                            "controls; nothing runs, no try is counted, no file"}
    miss = load(GUARD_FILE)
    if miss:
        return {"incomplete": True, "void": False, "missing": miss, "outcome": INCOMPLETE}
    g = body[GUARD_FILE]
    if g.get("pregate_sha256") != sha["pregate.json"]:
        why.append("guard.json was made under another pregate.json")
    common["guard_part"] = {k: g.get(k) for k in ("cells", "failed", "passes")}
    gmiss = [f for f in guard_files() if not (runs / f).exists()]
    if gmiss:
        return {"incomplete": True, "void": False, "missing": gmiss, "outcome": INCOMPLETE}
    gwhy: list[str] = []
    gsha, gprov = {}, {}
    grows = rows_of(runs, guard_files(), f"{GUARD_DIR}/", {c: (base_of(c),) for c in CELLS}, gwhy, gsha, gprov,
                    sha["pregate.json"])
    if gsha != g.get("files_sha256"):
        why.append("guard.json: its files' sha256 are not the guard files'")
    why += gwhy
    sha |= gsha
    prov |= gprov
    if g.get("passes") is not True:
        failed = [x for x in (g.get("failed") or []) if x in (BASE, BASE_X)]
        # row K26: "a void counts": the three tries, though the gate runs the bases only and stops before any try
        return {"incomplete": False, "void": True, **common, "files_sha256": sha, "provenance": prov,
                "tries_counted": VOID_TRIES, "stage0_counted": 1,
                "why": [f"the reproduction gate failed: {failed or 'unrecorded'}", *why],
                "outcome": "void: the reproduction gate failed before any try; no file; the three tries count"}
    missing = [f for f in score_files() if not (runs / f).exists()]
    if missing:
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    want = {c: (X_ARMS if c[0] == "x2" else G_ARMS) for c in CELLS}
    rows = rows_of(runs, score_files(), "", want, why, sha, prov, sha["pregate.json"])
    for f, p in prov.items():
        if p.get("tree_clean") is not True:
            why.append(f"{f}: a dirty tree")
    commits = {p.get("commit") for p in prov.values()}
    machines = {p.get("machine") for p in prov.values()}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(map(str, commits))}")
    if len(machines) > 1:
        why.append(f"results from different machines {sorted(map(str, machines))}")
    try:
        refs = guard_refs() if refs is None else refs
        if guard_eval(grows, refs)["passes"] is not True:
            why.append("the guard files, read again against the committed rows, do not pass the reproduction gate")
    except Refused as e:
        why.append(str(e))
    per: dict = {}
    match: dict = {}
    for c, rr in rows.items():
        if [q["seed"] for q in rr] != list(SEEDS):
            why.append(f"{c[0]} {c[1]}: not the ten seeds")
            continue
        arms = (harness.FLOOR, *want[c]) if c[0] != "dm" else want[c]
        per[c] = summarise(rr, arms, base_of(c))
        dd = k18.distance(rr, base_of(c), grows.get(c) or [], base_of(c))
        match[f"{c[0]} {c[1]}"] = dd
        if not k18.guard_passes(dd):
            why.append(f"{c[0]} {c[1]}: the scoring part's base is off the guard part's")
    dm_rows = {n: rows[("dm", n)] for n in TASKS}
    check = k18.dm_check(dm_rows) if all(len(v) == len(SEEDS) for v in dm_rows.values()) else None
    base = {"incomplete": False, **common, "files_sha256": sha, "provenance": prov, "guard_match": match,
            "dm_check_u3d_report_only": check, "tries_counted": MAX_TRIES, "stage0_counted": 1,
            "cells": {f"{c[0]} {c[1]}": v for c, v in per.items()}}
    if why:
        return base | {"void": True, "tries_counted": VOID_TRIES, "why": why,
                       "outcome": "void: no file; the three tries count; T3's built and listed files unchanged"}
    dec = decide(per)
    return base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                   "decision": dec, "outcome": dec["outcome"]}


def _pregate_report(pre: dict) -> dict:
    """The pre-gate's figures without gene-level values."""
    out = {"runs": pre.get("runs"), "dropped": pre.get("dropped"), "folds": {}}
    for f, v in (pre.get("folds") or {}).items():
        out["folds"][f] = {k: v.get(k) for k in ("train", "n_train", "rank", "lambda", "query", "rand_gene", "e",
                                                 "e_rand", "figures", "compared", "left_out", "corr", "dropped")}
    return out


def _fmt(v) -> str:
    return f"{v:.4g}" if isinstance(v, float) else str(v)


def text(res: dict) -> str:
    lines = [f"Try K26: {res['outcome']}."]
    if res.get("incomplete"):
        miss = res["missing"]
        return "\n".join([*lines, f"missing: {len(miss)} file(s), e.g. " + ", ".join(miss[:6])]) + "\n"
    lines.append(f"tries counted: {res['tries_counted']} (and stage 0, a counted diagnostic: {res['stage0_counted']})")
    s0 = res.get("stage0") or {}
    lines.append(f"stage 0: passes {s0.get('passes')}; held-out {s0.get('n_held')}; vs mean {s0.get('vs_mean')}; vs "
                 f"random queries {s0.get('vs_rand')}")
    for h, r in (s0.get("per_gene") or {}).items():
        lines.append(f"  {h}: rho map {_fmt(r.get('rho_map'))}, LOO mean {_fmt(r.get('rho_mean'))}, random "
                     f"{_fmt(r.get('rho_rand'))}")
    pre = res.get("pregate") or {}
    lines.append(f"pre-gate: runs {pre.get('runs')}; dropped {pre.get('dropped') or 'none'}")
    for f, v in (pre.get("folds") or {}).items():
        lines.append(f"  fold {f}: query {v.get('query')}, q_rand gene {v.get('rand_gene')}, trained on "
                     f"{v.get('train')} (rank {v.get('rank')}, lambda {v.get('lambda')}); e {_fmt(v.get('e'))}, "
                     f"e_rand {_fmt(v.get('e_rand'))}; corr {v.get('corr')}; figures {v.get('figures')}")
    lines.append(f"breadth of the queries' profiles (measured from T1): {res.get('breadth')}")
    lines.append(f"outside genes left out: {res.get('outside_left_out')}")
    gp = res.get("guard_part") or {}
    lines.append(f"reproduction gate: passes {gp.get('passes')}, failed {gp.get('failed') or 'none'}")
    for k, v in (gp.get("cells") or {}).items():
        lines.append(f"  {k}: passes {v['passes']}; {k21._dist(v.get('distance'))}")
    if res["void"]:
        lines.append("VOID: " + "; ".join(res["why"]))
    elif res.get("stopped"):
        lines.append(f"stopped at {res['stopped']}: no fold was read")
    else:
        dec = res["decision"]
        lines.append(f"commit {res['commit']}, machine {res['machine']}.")
        lines.append(f"pass line: (a) {dec['a']}, (b) {dec['b']}, (c) {dec['c']}, (d) {dec['d']}, (e) {dec['e']}; "
                     f"ahead in {dec['ahead_cells'] or 'none'}; behind in {dec['behind_cells'] or 'none'}; gated "
                     f"(arm, base) {dec['gated_arm_base']}")
    for k, d in (res.get("guard_match") or {}).items():
        lines.append(f"  scoring part's base against the guard part's, {k}: {k21._dist(d)}")
    if res.get("dm_check_u3d_report_only"):
        lines.append("U3d's check arm in the dm world, reported only: " + "; ".join(
            f"{n} {k21._dist(c['distance'])}" for n, c in res["dm_check_u3d_report_only"].items()))
    lines.append("reported, deciding nothing: the tables below.")
    for cell, v in (res.get("cells") or {}).items():
        lines.append(f"{cell}:")
        for a, s in v["arms"].items():
            lines.append(f"    {a:<9} {s['mean']:.2f} (sd {s['sd']:.2f}); gated {s['severity_gated']} "
                         f"{s['severity_states']}; skills " + ", ".join(f"{m} {x:.3f}" for m, x in s["skills"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
        for a, terms in v["per_term"].items():
            lines.append(f"    per-term skill, {a} - base: " + ", ".join(f"{m} {x:+.3f}" for m, x in terms.items())
                         + f"; de_score + de_direction {v['de_change'][a]:+.3f}")
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    if res.get("reading_code"):
        rc = res["reading_code"]
        lines.append(f"read by {rc['script']} {rc['script_sha256'][:12]} (context_map {rc['context_map_sha256'][:12]}) "
                     f"at commit {rc['commit']}, clean {rc['tree_clean']}, machine {rc['machine']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    head_ok(st, a.commit)
    out = Path(a.out)
    targets = [out / "reading.json", out / "reading.txt"]
    if any(p.exists() for p in targets):
        print(f"refused: {out} holds a reading already; it is never overwritten", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    try:
        res = reading(Path(a.runs))
        if not res.get("incomplete"):
            res["reading_code"] = reading_code()
            res["row"] = rid
        t = text(res)
        js = json.dumps(_jsonable(res), indent=1, default=float)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if res.get("incomplete"):
        print("INCOMPLETE: not a reading", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(js)
    (out / "reading.txt").write_text(t)
    print(f"wrote {out / 'reading.json'} and {out / 'reading.txt'}", flush=True)


def _jsonable(o):
    if isinstance(o, dict):
        return {(f"{k[0]} {k[1]}" if isinstance(k, tuple) else k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    return o


# ---------------------------------------------------------------------------------------------------------------------


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("features")
    s = sub.add_parser("stage0")
    p = sub.add_parser("pregate")
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("guard", "score"), default="score")
    r.add_argument("--world", choices=(*WORLDS, "x2"), required=True)
    r.add_argument("--task", choices=sorted(TASKS))
    r.add_argument("--seeds", default=None)
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--threads", type=int, default=8)
    for q in (f, s, p, r):
        q.add_argument("--out", required=True)
    for q in (f, p, r):
        q.add_argument("--data", default=str(harness.DATA_DIR))
    for q in (p, r):
        q.add_argument("--run-dir", default=str(k18.RUN_DIR))
        q.add_argument("--decision", default=str(k18.DECISION))
    gd = sub.add_parser("guard")
    gd.add_argument("--runs", required=True)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    for q in (f, s, p, r, gd, d):
        q.add_argument("--commit", required=True)
    a = ap.parse_args(argv)
    {"features": cmd_features, "stage0": cmd_stage0, "pregate": cmd_pregate, "run": cmd_run, "guard": cmd_guard,
     "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
