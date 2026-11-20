#!/usr/bin/env python
"""Try K23 (experiments/harness/TRIES.md; the row is placed before any run, and `pregate` and `run` refuse without it):
the barred factor times the same outside-knockout table as K22 (`all`, the same bytes and pins), raised to two larger
exponents, e = 2 and e = 4, read against the barred s1 on the two Gata4 replicates in K18's two practice worlds. A
labelled reading that may name at most one T3 candidate file by the rule fixed in the row; not the adoption rule. The
row follows a near miss on the same fold (K22), and says so.

Built on k22_outside_forms.py, outside_forms.py, k18_built_files.py and k21_factor_forms.py, all imported unchanged.
K22's table loader (`k22.load_ext`, which fills `k22._EXT`), its fold_factors (K18's plus the table in the world's
gene order and its permutation within deciles of the world's basis), its stamp check and its frozen-file check are
reused as they are. For a run, three of k18's module names are replaced inside this process only (`k18_arms`): ARMS,
build_arm and fold_factors (K22's), restored on exit and on an exception. So every arm draws the barred s1's cells
exactly as K18's barred_form did, and K18's barred_form is this try's base and first guard.

The frozen build applies a factor once, after the cells are drawn: t3_decompose.apply_gene_factor thins counts
binomially where the factor is below 1 and tops them up by Poisson where it is above 1; it neither renormalises nor
clips. build returns log1p CP10k. So the size twin's constant k_E changes no gene's share of a cell in expectation; it
acts through count depth only.

Arms (phi_M = transfer_mab21l2_all as loaded; T = the pooled table `all` as loaded; perm = outside_forms.permuted
within deciles of the basis K18's fold_factors receives in each world, seed 20261113, 10 bins):
  barred_form   phi_M                        base and guard 1 (K18's arm value for value; not counted inside it)
  x_e1          phi_M * T                    guard 2, the reproduction guard: K22's x_all through the same code path
                                             (not counted inside it; never named)
  x_e2, x_e4    phi_M * T ** E               the tries
  x_eE_perm     phi_M * perm(T) ** E         controls (perm(T) ** E = perm(T ** E))
  x_eE_size     k_E * phi_M ** a_E           size and depth twins (labelled controls, counted, never named), with
                                             a_E = SD(ln product_E) / SD(ln phi_M) and k_E = den_total(product_E) /
                                             den_total(phi_M ** a_E), computed once by the pre-gate (pregate.json)

  pregate --out DIR                                   the pre-gate: DIR/pregate.json and its lines (no score)
  run --world plain --task g1|g2 [--seeds 0-1] --out DIR   one plain task, seed by seed (resumable; reads
                                                      DIR/pregate.json)
  run --world dm [--seeds ...] [--workers 8 --threads 2] --out DIR   both replicates, row by row (resumable)
  read --runs DIR --out DIR                           the fixed reading: reading.json and reading.txt (never over an
                                                      existing one: exit 4)

Exit codes: 0 done; 1 a crash (masked traceback; seeds written so far stay); 2 refused (no row naming the script, a
dirty tree, pins not filled or not K22's, TABLES.sha256 untracked or off its pin, K22's committed pregate.json or
reading.json untracked, the folder variable unset or refused, bytes or records that disagree, the frozen T3 files off
their pins, a run without its pregate.json or with one that disagrees with the tables loaded or the world's basis); 3
INCOMPLETE (read: a file missing; nothing written); 4 read: the reading files exist already; 5 pregate: both
exponents dropped, nothing runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import tempfile
import traceback
from contextlib import contextmanager
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (the paired summary)
import k8_s2_form as k8  # noqa: E402  (the upload's central composition, for the record)
import k18_built_files as k18  # noqa: E402  (the worlds, tasks, stream, seeds, build path and scoring; unchanged)
import k21_factor_forms as k21  # noqa: E402  (the seed parser, the summary and K18's committed rows; unchanged)
import k22_outside_forms as k22  # noqa: E402  (the table loader, fold_factors, checks and constants; unchanged)
import outside_forms as ext  # noqa: E402  (the table and the product, shared with build_outside_form.py)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

SCRIPT = "experiments/t3_decomp/k23_outside_exponents.py"
ROW_SCRIPT = f"`{SCRIPT}`"  # the ledger row names this script
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"  # looked up at call time (row_present)
STREAM = k18.STREAM  # "t3_s2_form"
BASE = k18.BARRED  # "barred_form"
REPRO = "x_e1"
TABLE = "all"
EXPONENTS = (2.0, 4.0)
ARMS = (BASE, REPRO, "x_e2", "x_e2_perm", "x_e2_size", "x_e4", "x_e4_perm", "x_e4_size")
NAMEABLE = ("x_e2", "x_e4")
GUARDS = (BASE, REPRO)
CONTROL_OF = {"x_e2": "x_e2_perm", "x_e4": "x_e4_perm"}
SIZE_OF = {"x_e2": "x_e2_size", "x_e4": "x_e4_size"}
EXP_OF = {REPRO: 1.0, "x_e2": 2.0, "x_e2_perm": 2.0, "x_e2_size": 2.0, "x_e4": 4.0, "x_e4_perm": 4.0,
          "x_e4_size": 4.0}
TASKS = k18.TASKS
WORLDS = k18.WORLDS
MIN_GAIN, MIN_AHEAD, MAX_EXTRA_GATED = k22.MIN_GAIN, k22.MIN_AHEAD, k22.MAX_EXTRA_GATED  # 0.3, 8 of 10, +1
GUARD_MEAN, GUARD_SEED = k22.GUARD_MEAN, k22.GUARD_SEED  # 0.1 and 0.5 (K17's tolerance)
X_E1_TOL_MEAN = 0.1  # guard 2: x_e1's ten-seed mean against K22's committed x_all (K17's mean tolerance)
MAX_SD_LN_TABLE = 0.25  # (i) SD(ln T ** E) above this drops E
MAX_AT_CLIP = 45  # (ii) the product or its control in either world with more genes at or beyond abs ln 4 drops E
MAX_DEN_SHIFT = 0.05  # (iii) the product's den-weighted total off phi_M's by more than this (relative) drops E
MAX_SD_RATIO = 1.05  # (iv) SD(ln product) above this times SD(ln phi_M) drops E
DOSE_MIN = 0.0  # (f): x_eE - x_e1 (and x_e4 - x_e2) above this mean, with MIN_AHEAD seeds positive
PIVOT_BEHIND = -0.3  # naming: x_e4 - x_e2 below this mean in a cell keeps x_e2 (= -MIN_GAIN)
K23_LISTING_SHA256 = "884b294bb8ea15b4e97b8f9f3d72869c6f9a375c55a5ca11a6daf3c70b8d94ea"  # K22's, read by K22
K23_TABLE_SHA256 = "9490f524e40e146480c775f85ae8829ec9c005607104ad8ace129c5c46ed6228"  # `all`, read by K22
K18_TASKS = k22.K18_TASKS
K18_SHA256 = dict(k22.K18_SHA256)
RESULTS = HERE / "results" / "k23_outside_exponents"
K22_RESULTS = HERE / "results" / "k22_outside_forms"
REL_TOL = 1e-12  # pregate.json's exponent-1 figures against K22's committed pregate.json
EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_NOTHING_RUNS = 2, 3, 4, 5
_K23: dict = {}  # {E: {"a": a_E, "k": k_E}} from pregate.json, set before any arm builds
NOT_NAMED_E4 = "e = 4 alone qualified; not named: the smaller exponent does not point the same way"
NO_FILE = "no file: nothing qualifies; T3's list is unchanged"


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


def ekey(e: float) -> str:
    """The exponent as a pregate.json key: '2', '4'."""
    return f"{float(e):g}"


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row and the run's preconditions


def row_present(path: Path | None = None) -> str:
    """The id of the TRIES.md row that names this script; refuses (exit 2) when there is none or more than one."""
    path = TRIES_MD if path is None else path
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.startswith("| ") and ROW_SCRIPT in ln]
    ids = sorted({ln[2:].split(" ", 1)[0] for ln in rows})
    if len(ids) != 1:
        refuse(f"{'no' if not ids else 'more than one'} row of {Path(path).name} names {ROW_SCRIPT} ({ids}); the row "
               "is committed before any run")
    return ids[0]


DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing written. Row K23: commit or restore first; INCOMPLETE, not VOID")


def pins_are_k22s() -> None:
    """outside_forms' pins are K22's, the ones K23 reads (no new pin)."""
    if ext.K22_LISTING_SHA256 is None or ext.K22_TABLE_SHA256.get(TABLE) is None:
        refuse("the pins K22_LISTING_SHA256 and K22_TABLE_SHA256 are not filled (outside_forms.py): refused")
    if ext.K22_LISTING_SHA256 != K23_LISTING_SHA256 or ext.K22_TABLE_SHA256.get(TABLE) != K23_TABLE_SHA256:
        refuse("outside_forms.py's pins are not the listing and `all` sha256 K22 pinned and read: refused")


def preconditions() -> tuple[str, dict]:
    """Before anything loads: the row, a clean tree, the listing tracked and at its pin, the pins K22's, the frozen
    files, K22's committed pregate.json and reading.json tracked. Returns (row id, stamp)."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    try:
        ext.read_listing()
    except ext.TableRefused as e:
        refuse(str(e))
    pins_are_k22s()
    k22.frozen_ok()
    for f in ("pregate.json", "reading.json"):
        p = K22_RESULTS / f
        if not p.is_file() or not ext.tracked(p):
            refuse(f"K22's committed {f} is missing or not tracked by git: refused")
    return rid, st


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate


def _world_vectors(v: np.ndarray, genes: list[str] | None, world: dict) -> np.ndarray:
    """A vector in the barred factor's file order, put in a world's gene order (no reorder when no genes given)."""
    if genes is None or world.get("genes") is None:
        return np.asarray(v, dtype=np.float64)
    return k18.in_gene_order(v, list(genes), list(world["genes"]))


def pregate_exponents(phi_m, loaded: dict, bases_by_world: dict, den, genes: list[str] | None = None,
                      k22_pre: dict | None = None) -> dict:
    """The pre-gate of row K23, target-free, from phi_M, the loaded table, each world's decile basis and WT E8.75's
    mean within-cell proportions (den) alone; a pure function.

    bases_by_world: {world: {"genes": the world's gene order or None, "basis": that world's basis}}; phi_M, T and den
    are in the file order `genes`. k22_pre: outside_forms.pregate's output (computed here when None), kept whole under
    'k22'. An exponent is dropped (with its try, its control and its twin) when K22's pre-gate drops `all`, or at E:
    (i) SD(ln T ** E) above MAX_SD_LN_TABLE; (ii) the product, or its control in either world, with more than
    MAX_AT_CLIP genes at or beyond abs ln 4; (iii) the product's den-weighted total off phi_M's by more than
    MAX_DEN_SHIFT (relative); (iv) SD(ln product) above MAX_SD_RATIO times SD(ln phi_M). An E = 2 drop drops E = 4."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    pre22 = ext.pregate(phi_m, loaded, None, den) if k22_pre is None else k22_pre
    out: dict = {"k22": pre22, "n_genes": int(phi_m.size), "exponents": {}, "dropped": {},
                 "constants": {"exponents": list(EXPONENTS), "max_sd_ln_table": MAX_SD_LN_TABLE,
                               "max_at_clip": MAX_AT_CLIP, "max_den_shift": MAX_DEN_SHIFT,
                               "max_sd_ratio": MAX_SD_RATIO, "table": TABLE, "perm_seed": ext.PERM_SEED,
                               "perm_bins": ext.PERM_BINS, "worlds": sorted(bases_by_world)}}
    if TABLE in pre22["dropped"] or TABLE not in loaded["tables"]:
        why = f"K22's pre-gate drops `{TABLE}`: {pre22['dropped'].get(TABLE, 'not loaded')}"
        out["dropped"] = {ekey(e): why for e in EXPONENTS}
        return out
    T = np.asarray(loaded["tables"][TABLE][0], dtype=np.float64)
    sd_phi = float(np.std(np.log(phi_m)))
    den_phi = float(den @ phi_m / den.sum())
    earlier_dropped = False
    for e in EXPONENTS:
        Te = T ** e
        prod = ext.product(phi_m, 1.0, T, e)
        sd_prod = float(np.std(np.log(prod)))
        a = sd_prod / sd_phi
        k = float((den @ prod) / (den @ (phi_m ** a)))
        twin = k * phi_m ** a
        tf = ext.table_figures(phi_m, Te, den)
        pf = ext.figures(phi_m, prod, den)
        controls = {}
        for w in sorted(bases_by_world):
            wd = bases_by_world[w]
            phi_w = _world_vectors(phi_m, genes, wd)
            T_w = _world_vectors(T, genes, wd)
            den_w = _world_vectors(den, genes, wd)
            controls[w] = ext.figures(phi_w, ext.product(phi_w, 1.0, ext.permuted(T_w, wd["basis"]), e), den_w)
        why = []
        if tf["sd_ln"] > MAX_SD_LN_TABLE:
            why.append(f"SD(ln T ** {ekey(e)}) {tf['sd_ln']:.4f} (> {MAX_SD_LN_TABLE})")
        if pf["genes_at_clip"] > MAX_AT_CLIP:
            why.append(f"the product has {pf['genes_at_clip']} genes at or beyond the clip (> {MAX_AT_CLIP})")
        for w, c in controls.items():
            if c["genes_at_clip"] > MAX_AT_CLIP:
                why.append(f"the control in world {w} has {c['genes_at_clip']} genes at or beyond the clip "
                           f"(> {MAX_AT_CLIP})")
        shift = pf["den_total"] / den_phi - 1.0
        if abs(shift) > MAX_DEN_SHIFT:
            why.append(f"the product's den-weighted total moves {shift:+.2%} (beyond {MAX_DEN_SHIFT:.0%})")
        if sd_prod > MAX_SD_RATIO * sd_phi:
            why.append(f"SD(ln product) {sd_prod:.4f} exceeds {MAX_SD_RATIO} x SD(ln phi_M) {sd_phi:.4f}")
        if earlier_dropped:
            why.append("a smaller exponent is dropped")
        out["exponents"][ekey(e)] = {"exponent": float(e), "table": tf, "product": pf, "controls": controls,
                                     "size": ext.figures(phi_m, twin, den), "a": float(a), "k": float(k),
                                     "den_shift": float(shift), "sd_ln_phi_m": sd_phi, "sd_ln_product": sd_prod}
        if why:
            out["dropped"][ekey(e)] = "; ".join(why)
            earlier_dropped = True
    return out


def dropped_arms(pre: dict) -> dict[str, str]:
    """The pre-gate's exponent drops as arm drops; when every exponent drops, x_e1 is not run either."""
    out = {}
    for e in EXPONENTS:
        why = pre["dropped"].get(ekey(e))
        if why is not None:
            for a in (f"x_e{ekey(e)}", f"x_e{ekey(e)}_perm", f"x_e{ekey(e)}_size"):
                out[a] = f"exponent {ekey(e)}: {why}"
    if all(ekey(e) in pre["dropped"] for e in EXPONENTS):
        out[REPRO] = "every exponent dropped: nothing runs"
    return out


def active_arms(pre: dict) -> tuple[str, ...]:
    d = dropped_arms(pre)
    return tuple(a for a in ARMS if a not in d)


def nothing_runs(arms: tuple[str, ...]) -> bool:
    return not any(t in arms for t in NAMEABLE)


def constants_of(pre: dict) -> dict[str, dict]:
    """{E key: {"a", "k"}} of the exponents that passed (the twins' constants)."""
    return {k: {"a": v["a"], "k": v["k"]} for k, v in pre["exponents"].items() if k not in pre["dropped"]}


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def arm_factor(arm: str, F: dict) -> np.ndarray:
    """The factor an arm passes in place of phi, in the reference's gene order."""
    if arm == BASE:
        return F["phi_barred"]
    if arm not in EXP_OF:
        raise RuntimeError(f"unknown arm {arm!r}")
    e = EXP_OF[arm]
    if arm.endswith("_size"):
        c = _K23.get(ekey(e))
        if c is None:
            raise RuntimeError(f"the twin's constants for exponent {ekey(e)} are not loaded from pregate.json")
        return c["k"] * F["phi_barred"] ** c["a"]
    src = F["ext_perm"] if arm.endswith("_perm") else F["ext"]
    return ext.product(F["phi_barred"], 1.0, src[TABLE], e)


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """One arm's drawn cells, as K22's build_arm: every arm draws from stable_rng(seed, STREAM) through bs.construct
    with the barred spec's keywords; only the factor passed differs. barred_form passes the keywords unchanged."""
    if arm not in ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    if arm != BASE:
        kw["phi"] = arm_factor(arm, F)
    return bs.construct(ref, n, rng, **kw)


_K18_ORIGINAL = dict(k22._K18_ORIGINAL)  # k18's own ARMS, build_arm and fold_factors


@contextmanager
def k18_arms(arms: tuple[str, ...]):
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's (fold_factors is K22's),
    restored on exit and on an exception."""
    if not set(arms) <= set(ARMS) or BASE not in arms:
        raise RuntimeError(f"arms {arms} are not a set of this try's arms with the base")
    k18.ARMS, k18.build_arm, k18.fold_factors = tuple(arms), build_arm, k22.fold_factors
    k18._FOLD.clear()
    try:
        yield
    finally:
        k18.ARMS = _K18_ORIGINAL["ARMS"]
        k18.build_arm = _K18_ORIGINAL["build_arm"]
        k18.fold_factors = _K18_ORIGINAL["fold_factors"]
        k18._FOLD.clear()


def candidates(arms: tuple[str, ...]) -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: k18._build_plain(ctx, arm),
                                 params={"form": arm, "stream": STREAM, "composition": k8.spec()["composition"],
                                         "factor": ("the barred s1's" if arm == BASE else
                                                    f"k23 ({arm}; table {TABLE}, exponent {ekey(EXP_OF[arm])})"),
                                         "perm_seed": ext.PERM_SEED if arm.endswith("_perm") else None},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


# ---------------------------------------------------------------------------------------------------------------------
# Running


def progress_log(world: str, task: str):
    """The harness's per-seed line lists every arm's score; this prints the world, the task and the seed only."""
    import re

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"K23 {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


class _Mismatch(RuntimeError):
    pass


def load_pregate(out: Path, record: dict) -> tuple[dict, tuple[str, ...]]:
    """The run folder's pregate.json (written by `pregate` before any scoring part): its tables must be the loaded
    ones; fills _K23 with the twins' constants. Refuses (exit 2) when it is missing or disagrees."""
    p = Path(out) / "pregate.json"
    if not p.is_file():
        refuse(f"{p} is missing: the pre-gate runs first (`pregate --out {out}`)")
    body = json.loads(p.read_text())
    if body.get("ext") != record:
        refuse("pregate.json records other tables than the ones loaded: refused")
    pre = body["pregate"]
    arms = active_arms(pre)
    if list(arms) != body.get("arms"):
        refuse("pregate.json's arms are not its pre-gate's: refused")
    if nothing_runs(arms):
        refuse("the pre-gate drops both exponents: nothing runs")
    _K23.clear()
    _K23.update(constants_of(pre))
    return pre, arms


def _pregate_record(pre: dict, arms: tuple[str, ...]) -> dict:
    return {"dropped": dropped_arms(pre), "arms": list(arms), "ext": dict(k22._EXT["record"]),
            "exponents": constants_of(pre)}


def _world_figures(F: dict, arms: tuple[str, ...]) -> dict:
    """Per world, target-free: each arm's figures against phi_M in this world's order."""
    phi = F["phi_barred"]
    return {a: ext.figures(phi, arm_factor(a, F)) for a in arms if a != BASE}


def check_world(pre: dict, world: str, F: dict, arms: tuple[str, ...]) -> None:
    """The run's own controls equal pregate.json's controls in this world (the pre-gate's basis is the run's)."""
    for e in EXPONENTS:
        a = f"x_e{ekey(e)}_perm"
        if a not in arms:
            continue
        want = pre["exponents"][ekey(e)]["controls"].get(world)
        got = ext.figures(F["phi_barred"], arm_factor(a, F))
        if want is None or any(not _close(got[k], want.get(k)) for k in got):
            raise _Mismatch(f"the {world} world's control at exponent {ekey(e)} is not the pre-gate's (the pre-gate's "
                            "basis is not this world's): no task file written")


def _close(a, b) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        return (a is not None and b is not None and math.isclose(float(a), float(b), rel_tol=REL_TOL, abs_tol=0.0)) \
            or a == b
    return a == b


def run_plain(a, st: dict, out: Path) -> None:
    folder = out / f"plain_{a.task}"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    k22.load_ext(k18._BARRED["genes"], k18._BARRED["phi"])
    pre, arms = load_pregate(out, k22._EXT["record"])
    with k18_arms(arms):
        for s in k21.parse_seeds(a.seeds):
            path = folder / f"seed_{s}.json"
            if path.exists():
                print(f"K23 plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            res = harness.run_task(TASKS[a.task], candidates(arms), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log("plain", a.task))
            got = res["provenance"]
            if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
                raise k22._Unclean(f"the harness's stamp says commit {got.get('commit')}, clean "
                                   f"{got.get('tree_clean')}: no task file written (INCOMPLETE)")
            k22._stamp_ok(st)
            fold = k18._FOLD.get("F")
            if fold is None:
                raise _Mismatch("the fold's factors were not made: no task file written")
            check_world(pre, "plain", fold, arms)
            res["provenance"]["machine"] = st["machine"]
            res |= {"world": "plain", "task_source": SCRIPT, "seed": s, "pregate": _pregate_record(pre, arms),
                    "factor_record": k7._PREP["fold"]["record"], "barred_factor": rec,
                    "world_figures": _world_figures(fold, arms)}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)


def run_dm(a, st: dict, out: Path) -> None:
    """The depth-matched world, set up exactly as k22.run_dm sets it up, then each (replicate, seed) row through
    k18.dm_row with this try's arms; rows already written at this commit on this machine are skipped."""
    import multiprocessing as mp

    import orth_diag as od

    folder = out / "dm"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    seeds = k21.parse_seeds(a.seeds)
    bs = harness._t3_head()
    t3d = bs.t3d
    spec, dec, rules = k18.barred_spec(Path(a.decision))
    if list(rules["protocol"]["seeds"]) != list(SEEDS):
        raise RuntimeError("the frozen script's seeds are not vec.adoption.SEEDS")
    paths = od.context_paths(Path(a.data), with_targets=True)
    n_par = a.workers * a.threads
    by_name = {v: k for k, v in k18.DM_TARGETS.items()}
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=Path(tmp))
        if sorted(T.name for T in ctx.targets) != sorted(k18.DM_TARGETS.values()):
            raise RuntimeError(f"targets {[T.name for T in ctx.targets]}, expected {sorted(k18.DM_TARGETS.values())}")
        k18.quiet_selfcheck(t3d, ctx)
        genes = list(ctx.genes)
        phi_b, frec = k18.barred_factor(Path(a.run_dir), dec, spec, genes)
        k22.load_ext(genes, phi_b)
        pre, arms = load_pregate(out, k22._EXT["record"])
        ref_C, ref_lab, ref_P, gap = od.reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        wt_path = Path(paths["wt"])
        series = bs.wt_revert_series(wt_path, {n: Path(a.data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                               if n != wt_path.stem})
        phi, per, group, rrec = bs.wt_revert_factors(wt_path, series, genes, ctx.wt.pool.C)
        refw = bs.Reference(preset=ctx.preset_name, rules=rules, genes=genes, V=ctx.V, pool=ctx.wt.pool,
                            coords=np.zeros((ctx.wt.pool.C.shape[0], 3)))
        files = [*(Path(a.data) / f for f in od.WT_FILES.values()), *(Path(a.data) / f for f in od.KO_FILES.values()),
                 *series.values(), Path(a.decision), Path(a.run_dir) / "factors" / f"{od.BARRED}.csv"]
        setup = {"world": "dm", "barred_factor": frec, "wt_revert": {k: v for k, v in rrec.items() if k != "series"},
                 "reference_counts_largest_gap": gap, "pregate": _pregate_record(pre, arms),
                 "files_sha256": {str(f): t3d.sha256(f) for f in dict.fromkeys(files)}}
        sp = folder / "setup.json"
        if sp.exists():
            old = json.loads(sp.read_text())
            if {k: old.get(k) for k in setup} != json.loads(json.dumps(setup, default=float)):
                sys.exit(f"{sp} records another set-up (inputs, factor, tables or pre-gate): not run")
        else:
            k18.write_task_file(sp, setup | {"provenance": st})
        with k18_arms(arms):
            F = k18.fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b,
                                 spec)
            check_world(pre, "dm", F, arms)
            wf = _world_figures(F, arms)
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder / f"{by_name[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                k22._stamp_ok(st)
                path = folder / f"{by_name[row['target']]}_seed_{row['seed']}.json"
                k18.write_task_file(path, row | {"world": "dm", "task": by_name[row["target"]], "arms": list(arms),
                                                 "pregate": _pregate_record(pre, arms), "provenance": st,
                                                 "world_figures": wf, "task_source": SCRIPT})
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
    _, st = preconditions()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        (run_plain if a.world == "plain" else run_dm)(a, st, out)
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seeds written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)


def _fmt(v) -> str:
    return f"{v:.6g}" if isinstance(v, float) else str(v)


def _figs(r: dict) -> str:
    return ", ".join(f"{k} {_fmt(v)}" for k, v in r.items())


def pregate_lines(pre: dict) -> list[str]:
    lines = ["pre-gate (the table, phi_M and the worlds' bases alone; no score):"]
    lines += ["  K22's " + ln.lstrip() for ln in k22.pregate_lines(pre["k22"])]
    for k, v in pre["exponents"].items():
        lines.append(f"  e = {k}: a_E {v['a']!r}, k_E {v['k']!r}, den shift {v['den_shift']:+.4%}")
        lines.append(f"    table ** {k}: {_figs(v['table'])}")
        lines.append(f"    product: {_figs(v['product'])}")
        for w, c in v["controls"].items():
            lines.append(f"    control ({w} basis): {_figs(c)}")
        lines.append(f"    size twin: {_figs(v['size'])}")
    lines.append("  exponents dropped: " + ("; ".join(f"e = {k} ({v})" for k, v in pre["dropped"].items())
                                            or "none"))
    return lines


def world_bases(data: Path, genes: list[str], bs) -> dict:
    """Each world's decile basis, as K18's fold_factors receives it, from WT E8.75 alone (no target file is read):
    plain, harness.t3_mean_expression of the pool the harness loads (the panel board's genes, duplicates dropped);
    dm, the mean of the frozen script's WT log1p CP10k (build_context's Xw), in its gene order."""
    from vec.prep import load_stage

    boards = harness.load_boards()
    tasks = list(TASKS.values())
    keys = {(t.panel_board, t.reference, harness.base_role(t)) for t in tasks}
    if len(keys) != 1 or any(t.target_filter is not None for t in tasks):
        raise RuntimeError("g1 and g2 do not share one reference pool and panel: the plain basis is not one")
    (panel, ref_file, _), = keys
    pool = load_stage(Path(data) / ref_file, genes=boards[panel].genes, drop_duplicates=True)
    plain = {"genes": list(pool.var_names), "basis": harness.t3_mean_expression(pool)}
    wt = load_stage(Path(data) / "T2" / "E8.75.h5ad", drop_duplicates=True)
    Xw = bs.t3d.lognorm(bs.t3d.dense_counts(wt)).astype(np.float64)
    dm = {"genes": list(wt.var_names), "basis": Xw.mean(0)}
    for w, d in (("plain", plain), ("dm", dm)):
        if sorted(d["genes"]) != sorted(genes):
            raise RuntimeError(f"the {w} world's genes are not the barred factor's")
    return {"plain": plain, "dm": dm}


def cmd_pregate(a) -> None:
    """The pre-gate before any scoring: phi_M through the builder's loader (WT E8.75's gene order), the table through
    outside_forms (K22's loader), K22's figures as K22 made them, each world's basis from WT E8.75 and the den-weighted
    totals on its mean within-cell proportions; DIR/pregate.json."""
    rid, st = preconditions()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "pregate.json"
    if path.exists():
        refuse(f"{path} exists")
    try:
        from vec.prep import load_stage

        bs = harness._t3_head()
        rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
        phi_m, genes = k18._BARRED["phi"], k18._BARRED["genes"]
        wt = load_stage(Path(a.data) / "T2" / "E8.75.h5ad", drop_duplicates=True)
        if list(wt.var_names) != list(genes):
            raise RuntimeError("WT E8.75's genes are not the barred factor's")
        C = np.asarray(bs.t3d.dense_counts(wt), dtype=np.float64)
        den = (C / np.maximum(C.sum(1, keepdims=True), 1)).mean(0)
        pre22 = k22.load_ext(genes, phi_m, harness.t3_mean_expression(wt), den)
        bases = world_bases(Path(a.data), list(genes), bs)
        pre = pregate_exponents(phi_m, k22._EXT["loaded"], bases, den, list(genes), pre22)
        arms = active_arms(pre)
        res = {"pregate": pre, "arms": list(arms), "dropped_arms": dropped_arms(pre), "row": rid,
               "barred_factor": rec, "ext": k22._EXT["record"], "table_records": k22._EXT["table_records"],
               "provenance": st}
        k18.write_task_file(path, res)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print("\n".join(pregate_lines(pre)), flush=True)
    if nothing_runs(arms):
        print("both exponents dropped: nothing runs and nothing is counted", flush=True)
        sys.exit(EXIT_NOTHING_RUNS)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _cells():
    return [(w, n) for w in WORLDS for n in TASKS]


def extra_pairs(arms: tuple[str, ...]) -> list[tuple[str, str]]:
    """The pairs beyond K22's (each arm against barred_form, each try against its permuted control): each try against
    its size twin and against x_e1, and x_e4 against x_e2."""
    out = [(t, SIZE_OF[t]) for t in NAMEABLE if t in arms and SIZE_OF[t] in arms]
    out += [(t, REPRO) for t in NAMEABLE if t in arms and REPRO in arms]
    if "x_e4" in arms and "x_e2" in arms:
        out.append(("x_e4", "x_e2"))
    return out


def _summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """k21.summarise with this try's control pairs, plus extra_pairs (k7._paired, as k21 pairs)."""
    saved = dict(k21.CONTROL_OF)
    k21.CONTROL_OF.clear()
    k21.CONTROL_OF.update(CONTROL_OF)
    try:
        out = k21.summarise(rows, arms)
    finally:
        k21.CONTROL_OF.clear()
        k21.CONTROL_OF.update(saved)
    for x, y in extra_pairs(arms):
        out["paired"][f"{x} - {y}"] = k7._paired(rows, x, y)
    return out


def _beats(v: dict | None) -> bool:
    return bool(v is not None and v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD)


def _dose(v: dict | None) -> bool:
    return bool(v is not None and v["mean"] > DOSE_MIN and v["ahead"] >= MIN_AHEAD)


def _status(per: dict, t: str, arms: tuple[str, ...]) -> dict:
    """The qualification (a) to (f) of one try."""
    if t not in arms:
        return {"ran": False, "qualified": False}
    d = {(w, n): per[w][n]["paired"][f"{t} - {BASE}"] for w, n in _cells()}
    ahead = {k: _beats(v) for k, v in d.items()}
    behind = {k: bool(v["mean"] < -MIN_GAIN) for k, v in d.items()}
    gated = {k: (per[k[0]][k[1]]["arms"][t]["severity_gated"], per[k[0]][k[1]]["arms"][BASE]["severity_gated"])
             for k in d}
    pair = lambda k, other: per[k[0]][k[1]]["paired"].get(f"{t} - {other}")  # noqa: E731
    a = (any(ahead[("plain", n)] for n in TASKS) and any(ahead[("dm", n)] for n in TASKS)
         and any(ahead[(w, "g1")] for w in WORLDS) and any(ahead[(w, "g2")] for w in WORLDS)
         and any(ahead[("dm", n)] and gated[("dm", n)][0] == gated[("dm", n)][1] for n in TASKS))
    b = not any(behind.values())
    cc = all(g <= g0 + MAX_EXTRA_GATED for g, g0 in gated.values())
    up = [k for k in d if ahead[k]]
    dd = (all(_beats(pair(k, CONTROL_OF[t])) for k in up)
          and any(k[0] == "dm" and _beats(pair(k, CONTROL_OF[t])) for k in up))
    ee = all(_beats(pair(k, SIZE_OF[t])) for k in up)
    ff = all(_dose(pair(k, REPRO)) for k in up)
    if t == "x_e4":
        ff = ff and all(_dose(pair(k, "x_e2")) for k in up)
    return {"ran": True, "qualified": bool(a and b and cc and dd and ee and ff), "a": bool(a), "b": bool(b),
            "c": bool(cc), "d": bool(dd), "e": bool(ee), "f": bool(ff),
            "ahead_cells": [f"{w} {n}" for (w, n), v in ahead.items() if v],
            "behind_cells": [f"{w} {n}" for (w, n), v in behind.items() if v],
            "means": {f"{w} {n}": float(v["mean"]) for (w, n), v in d.items()},
            "mean_of_four": float(np.mean([v["mean"] for v in d.values()])),
            "gated_vs_base": {f"{w} {n}": list(v) for (w, n), v in gated.items()}}


def decide(per: dict, arms: tuple[str, ...]) -> dict:
    """The rule fixed in the row, on the summaries per (world, task): at most one name."""
    status = {t: _status(per, t, arms) for t in NAMEABLE}
    q = [t for t in NAMEABLE if status[t]["qualified"]]
    pivot = None
    if q == ["x_e2", "x_e4"]:
        p = {(w, n): per[w][n]["paired"]["x_e4 - x_e2"] for w, n in _cells()}
        up = {k: _beats(v) for k, v in p.items()}
        pivot = {"ahead_cells": [f"{w} {n}" for (w, n), v in up.items() if v],
                 "below_cells": [f"{w} {n}" for (w, n), v in p.items() if v["mean"] < PIVOT_BEHIND],
                 "means": {f"{w} {n}": float(v["mean"]) for (w, n), v in p.items()}}
        pivot["e4"] = bool(any(up[("plain", n)] for n in TASKS) and any(up[("dm", n)] for n in TASKS)
                           and not pivot["below_cells"])
        name, rule = ("x_e4", "both qualified; x_e4 - x_e2 meets the bar") if pivot["e4"] else \
            ("x_e2", "both qualified; x_e4 - x_e2 does not meet the bar")
    elif q == ["x_e2"]:
        name, rule = "x_e2", "only x_e2 qualified"
    elif q == ["x_e4"]:
        s2 = status["x_e2"]
        consistent = bool(
            s2["ran"] and not s2["behind_cells"]
            and all(g <= g0 + MAX_EXTRA_GATED for g, g0 in s2["gated_vs_base"].values())
            and all(s2["means"][c] > 0 for c in status["x_e4"]["ahead_cells"]))
        name, rule = ("x_e4", "only x_e4 qualified; x_e2 points the same way") if consistent else \
            (None, NOT_NAMED_E4)
    else:
        name, rule = None, "nothing qualifies"
    names = [name] if name else []
    e = None if name is None else ekey(EXP_OF[name])
    outcome = (f"names {name}: phi_M ** x * T ** {e} on base B (build_outside_form.py --form xall-<B> "
               f"--table-exponent {float(e):.1f}), after (q)'s C list and after any K21 or K22 file, under the "
               "R14/R17 conditions" if name else NOT_NAMED_E4 if rule == NOT_NAMED_E4 else NO_FILE)
    return {"status": status, "qualified": q, "pivot": pivot, "names": names, "rule": rule,
            "files": [f"xall-B --table-exponent {float(e):.1f} (base B of decision log 2026-11-20 (q))"] if name
            else [], "outcome": outcome}


def expected_files() -> list[str]:
    return k21.expected_files()


INCOMPLETE = ("not a reading: pregate.json or a seed's task file is missing (a crash, a refusal such as a dirty tree, "
              "out of memory or a kill); nothing was read. Row K23: if the cause is outside the code, rerun the "
              "missing seeds into this folder at the run's commit on the same machine (the runs skip what is there); "
              "after a committed code fix, rerun everything at the new commit into a fresh folder and never read this "
              "one")


def pinned_ext() -> dict:
    return {"listing_sha256": K23_LISTING_SHA256, "tables": {TABLE: K23_TABLE_SHA256}}


def reading_code() -> dict:
    """Which code made the reading: this script's, outside_forms', K22's and K21's sha256, the commit, clean flag."""
    st = k18.stamp()
    h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
    return {"script": SCRIPT, "script_sha256": h(Path(__file__)), "outside_forms_sha256": h(HERE / "outside_forms.py"),
            "k22_sha256": h(HERE / "k22_outside_forms.py"), "k21_sha256": h(HERE / "k21_factor_forms.py"),
            "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}


def _same_figures(got: dict | None, want: dict | None) -> bool:
    if not isinstance(got, dict) or not isinstance(want, dict) or set(got) != set(want):
        return False
    return all(_close(got[k], want[k]) for k in want)


def k22_record(k22_dir: Path) -> tuple[dict | None, dict | None, list[str]]:
    """K22's committed pregate.json and reading.json (never its run folders)."""
    why = []
    out = []
    for f in ("pregate.json", "reading.json"):
        p = Path(k22_dir) / f
        if not p.is_file():
            why.append(f"K22's committed {f} is missing")
            out.append(None)
            continue
        out.append(json.loads(p.read_text()))
    return out[0], out[1], why


def guard2(rows: list[dict], k22_cell: dict | None, gated: int) -> dict:
    """x_e1's ten-seed mean and gated count against K22's committed x_all in one cell."""
    m = float(np.mean([r["scores"][REPRO] for r in rows]))
    if not k22_cell:
        return {"passes": False, "mean": m, "k22_mean": None, "gated": gated, "k22_gated": None}
    d = m - float(k22_cell["mean"])
    return {"passes": bool(abs(d) <= X_E1_TOL_MEAN and gated == k22_cell["severity_gated"]), "mean": m,
            "k22_mean": float(k22_cell["mean"]), "mean_diff": d, "gated": gated,
            "k22_gated": k22_cell["severity_gated"]}


def reading(runs: Path, refs: Path = K18_TASKS, u3d: Path = k18.U3D_RECORDS, k22_dir: Path = K22_RESULTS) -> dict:
    """Apply the reading of row K23 to a run folder."""
    runs = Path(runs)
    missing = [f for f in expected_files() if not (runs / f).exists()]
    if missing:
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    why: list[str] = []
    sha, prov = {}, {}
    pre_file, sha["pregate.json"] = k21._read_json(runs / "pregate.json")
    prov["pregate.json"] = pre_file["provenance"]
    pre = pre_file["pregate"]
    arms = active_arms(pre)
    if list(arms) != pre_file["arms"]:
        why.append("pregate.json: its arms are not the pre-gate's")
    pins = pinned_ext()
    erec = pre_file.get("ext") or {}
    if {"listing_sha256": erec.get("listing_sha256"), "tables": erec.get("tables")} != pins:
        why.append("pregate.json: its tables are not the pinned ones (listing or table sha256)")
    k22_pre, k22_read, kwhy = k22_record(k22_dir)
    why += kwhy
    if k22_read is not None:
        ke = k22_read.get("ext") or {}
        if {"listing_sha256": ke.get("listing_sha256"), "tables": ke.get("tables")} != pins or \
                ke.get("combine_sha256") != erec.get("combine_sha256"):
            why.append("pregate.json: its tables are not the ones K22's committed reading.json records")
    if k22_pre is not None:
        mine, theirs = pre.get("k22") or {}, k22_pre.get("pregate") or {}
        for kind in ("tables", "products"):
            if not _same_figures((mine.get(kind) or {}).get(TABLE), (theirs.get(kind) or {}).get(TABLE)):
                why.append(f"pregate.json: K22's {kind[:-1]} figures for `{TABLE}` differ from K22's committed "
                           f"pregate.json (relative {REL_TOL})")
    want_rec = {"dropped": dropped_arms(pre), "arms": list(arms), "ext": erec, "exponents": constants_of(pre)}
    rows: dict = {(w, n): [] for w, n in _cells()}
    for f in expected_files()[1:]:
        r, sha[f] = k21._read_json(runs / f)
        prov[f] = r["provenance"]
        world, rest = ("plain", f.split("/")[0][len("plain_"):]) if f.startswith("plain_") else ("dm", f[3:5])
        seed = int(f.rsplit("_", 1)[1][:-len(".json")])
        rec = r.get("pregate") or {}
        if rec.get("arms") != want_rec["arms"] or rec.get("dropped") != want_rec["dropped"]:
            why.append(f"{f}: arms or drops differ from pregate.json's")
        if rec.get("exponents") != want_rec["exponents"]:
            why.append(f"{f}: exponents, a_E or k_E differ from pregate.json's")
        if rec.get("ext") != erec:
            why.append(f"{f}: its tables differ from pregate.json's")
        if world == "plain":
            if r["task"]["key"] != TASKS[rest].key or [q["seed"] for q in r["rows"]] != [seed]:
                why.append(f"{f}: another task or seed")
                continue
            row = r["rows"][0]
            want = (harness.FLOOR, *arms)
        else:
            if r.get("task") != rest or r.get("seed") != seed:
                why.append(f"{f}: another task or seed")
                continue
            row = r
            want = arms
        miss = [x for x in want if x not in row.get("scores", {}) or x not in row.get("skills", {})
                or x not in row.get("t3", {})]
        if miss:
            why.append(f"{f}: missing arms {miss}")
            continue
        rows[(world, rest)].append(row)
    for f, p in prov.items():
        if p.get("tree_clean") is not True:
            why.append(f"{f}: a dirty tree")
    commits = {p.get("commit") for p in prov.values()}
    machines = {p.get("machine") for p in prov.values()}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(map(str, commits))}")
    if len(machines) > 1:
        why.append(f"results from different machines {sorted(map(str, machines))}")
    saved = dict(k21.K18_SHA256)
    k21.K18_SHA256.clear()
    k21.K18_SHA256.update(K18_SHA256)
    try:
        ref_rows, kwhy = k21.k18_rows(refs)
    finally:
        k21.K18_SHA256.clear()
        k21.K18_SHA256.update(saved)
    why += kwhy
    per: dict = {"plain": {}, "dm": {}}
    guard: dict = {}
    g2: dict = {}
    for (w, n), rr in rows.items():
        rr = sorted(rr, key=lambda q: q["seed"])
        if [q["seed"] for q in rr] != list(SEEDS):
            why.append(f"{w} {n}: not the ten seeds")
            continue
        per[w][n] = _summarise(rr, (harness.FLOOR, *arms) if w == "plain" else arms)
        if (w, n) in ref_rows:
            d = k18.distance(rr, BASE, ref_rows[(w, n)], BASE)
            guard[f"{w} {n}"] = d
            if not k18.guard_passes(d):
                why.append(f"guard 1: {w} {n}: barred_form is off K18's committed barred_form" + (
                    f" by {d['mean_diff']:+.3f} in the mean and {d['max_abs_diff']:.3f} in one seed" if d else
                    " (a seed is missing)") + f" (limits {GUARD_MEAN} and {GUARD_SEED})")
        if REPRO in arms:
            cell = (((k22_read or {}).get("worlds") or {}).get(w) or {}).get(n) or {}
            g = guard2(rr, (cell.get("arms") or {}).get("x_all"), per[w][n]["arms"][REPRO]["severity_gated"])
            g2[f"{w} {n}"] = g
            if not g["passes"]:
                why.append(f"guard 2: {w} {n}: x_e1 is off K22's committed x_all (mean {g['mean']:.2f} against "
                           f"{g['k22_mean']}, gated {g['gated']} against {g['k22_gated']}; limit {X_E1_TOL_MEAN})")
    check = k18.dm_check({n: sorted(rows[("dm", n)], key=lambda q: q["seed"]) for n in TASKS}, u3d)
    g1_ok = len(guard) == 4 and all(k18.guard_passes(d) for d in guard.values())
    g2_ok = REPRO not in arms or (len(g2) == 4 and all(g["passes"] for g in g2.values()))
    tries = sum(a not in GUARDS for a in arms) + (0 if g1_ok else 1) + (0 if g2_ok or REPRO not in arms else 1)
    base = {"incomplete": False, "files_sha256": sha, "provenance": prov, "arms": list(arms), "pregate": pre,
            "ext": erec, "guard": guard, "guard2": g2, "dm_check_u3d_report_only": check, "tries_counted": tries}
    if why:
        return base | {"void": True, "why": why, "worlds": per, "names": [],
                       "outcome": "void: no file; T3's list is unchanged"}
    return (base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                    "worlds": per} | decide(per, arms))


def text(res: dict) -> str:
    lines = [f"Try K23: {res['outcome']}."]
    if res.get("incomplete"):
        miss = res["missing"]
        return "\n".join([*lines, f"missing: {len(miss)} file(s), e.g. " + ", ".join(miss[:6])]) + "\n"
    lines.append(f"arms run: {', '.join(res['arms'])}; tries counted: {res['tries_counted']}")
    lines.append(f"tables: listing {res['ext'].get('listing_sha256')}, " + ", ".join(
        f"{n} {s}" for n, s in (res["ext"].get("tables") or {}).items()))
    lines += pregate_lines(res["pregate"])
    if res["void"]:
        lines.append("VOID: " + "; ".join(res["why"]))
    else:
        lines.append(f"commit {res['commit']}, machine {res['machine']}.")
        for t, s in res["status"].items():
            if not s["ran"]:
                lines.append(f"  {t}: not run (dropped by the pre-gate)")
                continue
            lines.append(f"  {t}: qualified {'yes' if s['qualified'] else 'no'} ((a) {s['a']}, (b) {s['b']}, (c) "
                         f"{s['c']}, (d) {s['d']}, (e) {s['e']}, (f) {s['f']}); ahead in {s['ahead_cells'] or 'none'}"
                         f"; behind in {s['behind_cells'] or 'none'}; mean of the four paired means "
                         f"{s['mean_of_four']:+.2f}; gated (arm, base) {s['gated_vs_base']}")
        lines.append(f"  qualified: {res['qualified'] or 'none'}; naming rule: {res['rule']}; pivot x_e4 - x_e2: "
                     f"{res['pivot'] if res['pivot'] is not None else 'not applied'}; named: {res['files'] or 'none'}")
    lines.append("guard 1 (barred_form against K18's committed barred_form): " + "; ".join(
        f"{k} {k21._dist(d)}" for k, d in res["guard"].items()))
    lines.append("guard 2 (x_e1 against K22's committed x_all): " + ("; ".join(
        f"{k} mean {g['mean']:.2f} against {g['k22_mean']}, gated {g['gated']} against {g['k22_gated']}"
        for k, g in (res.get("guard2") or {}).items()) or "not run"))
    lines.append("U3d's check arm in the dm world, reported only: " + "; ".join(
        f"{n} {k21._dist(c['distance'])}" for n, c in (res.get("dm_check_u3d_report_only") or {}).items()))
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    lines.append("g1 and g2 are two replicates of one knockout run; the controls, the twins and x_e1 never qualify.")
    for w in WORLDS:
        lines.append(f"\nWORLD {w} ({k18.WORLD_NAMES[w]}):")
        for n, t in TASKS.items():
            v = res["worlds"].get(w, {}).get(n)
            if v is None:
                lines.append(f"{w} {n} {t.key}: not read")
                continue
            lines.append(f"{w} {n} {t.key}:")
            for a, s in v["arms"].items():
                lines.append(f"    {a:<12} {s['mean']:.2f} (sd {s['sd']:.2f}); severity scored {s['severity_scored']}"
                             f", gated {s['severity_gated']} {s['severity_states']}, median beta "
                             f"{k21._f(s['median_beta'])}, r2 {k21._f(s['median_r2'])}; skills "
                             + ", ".join(f"{m} {x:.3f}" for m, x in s["skills"].items()))
            for k, d in v["paired"].items():
                lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
            for a, terms in v["per_term"].items():
                lines.append(f"    per-term skill, {a} - {BASE}: " + ", ".join(f"{m} {x:+.3f}" for m, x in
                                                                              terms.items()))
    if res.get("reading_code"):
        rc = res["reading_code"]
        lines.append(f"\nread by {rc['script']} {rc['script_sha256'][:12]} (outside_forms "
                     f"{rc['outside_forms_sha256'][:12]}, k22 {rc['k22_sha256'][:12]}, k21 {rc['k21_sha256'][:12]}) "
                     f"at commit {rc['commit']}, clean {rc['tree_clean']}, machine {rc['machine']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    targets = [out / "reading.json", out / "reading.txt"]
    if any(p.exists() for p in targets):
        print(f"refused: {out} holds a reading already; it is never overwritten (move it, or read into a fresh "
              "folder)", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    try:
        res = reading(Path(a.runs))
        if not res.get("incomplete"):
            res["reading_code"] = reading_code()
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
    except Exception:  # noqa: BLE001  (masked traceback, nothing written)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if res.get("incomplete"):
        print(t, end="")
        sys.exit(EXIT_INCOMPLETE)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(js)
    (out / "reading.txt").write_text(t)
    print(t, end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pregate")
    p.add_argument("--out", required=True)
    r = sub.add_parser("run")
    r.add_argument("--world", choices=WORLDS, required=True)
    r.add_argument("--task", choices=sorted(TASKS))
    r.add_argument("--seeds", default=None)
    r.add_argument("--out", required=True)
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--threads", type=int, default=8)
    for q in (p, r):
        q.add_argument("--data", default=str(harness.DATA_DIR))
        q.add_argument("--run-dir", default=str(k18.RUN_DIR))
        q.add_argument("--decision", default=str(k18.DECISION))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"pregate": cmd_pregate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
