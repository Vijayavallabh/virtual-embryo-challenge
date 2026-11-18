#!/usr/bin/env python
"""Try K18 (experiments/harness/TRIES.md, registered before any run): T3's two built last-upload files read against
each other on the two Gata4 replicates, in two practice worlds. A labelled reading for a ruling, not the adoption rule.
Whatever it shows, the try changes no entry; at most a priced question goes to the captain, who alone can change
T3's last upload.

Why. T3's last upload is the refused branch's s1 with the wild-type-only factor (tries K7b and K8b; docs/S2_RULES.md,
"T3: the second upload also carries a wild-type-only factor"). The reference list also holds the barred branch's s1:
the Mab21l2 factor (measured against WT E8.75) at exponent 1 with the same central composition and no depth step
(docs/S2_RULES.md, "Branch names"). On the Gata4 replicates the barred factor without the composition read 59.40 and
62.15 (gata4_records.csv) and the refused form 57.14 and 55.14 (K8b), but the barred file's own recipe, factor plus
composition, was never scored. In the depth-matched world of orth_diag.py (try U3d: the reference and the wild-type
pool thinned to each replicate's per-label depth) a stand-in for the first upload had its severity gated in every row,
as at the test (it scored 53.89 there, not the test's figure), so the two files are read there as well as in K8's plain
world.

Arms, the same in both worlds and all from one stream per world (STREAM, K8's), so every arm of the upload's form holds
the same drawn wild-type cells; every arm goes through the builder's own functions (build_submission.composition,
construct, prediction_inputs, load_factor, and on the depth-matched world wt_revert_factors):
  refused_form     the last upload as built: the central composition, then each drawn cell moved by its lineage
                   group's wild-type-only factor at exponent 0.5 (K8b's arm; on the plain world value for value)
  barred_form      the barred s1's recipe exactly as build_submission.py gives it for --branch barred --submission s1:
                   the spec by resolve_spec, the factor file by load_factor (sha256-checked against the decision file),
                   construct with prediction_inputs' keywords (the Mab21l2 factor at exponent 1, the central
                   composition, no depth ratio, no stage mixing)
  Labelled diagnostics, deciding nothing (read as paired score and median r2 against refused_form):
  resid_revert_k7  K7's form (the reference's own composition): each cell moved by its group's wild-type-only factor
                   with its component along the pooled factor removed (residual_factors), exponent 0.5
  resid_perm_k7    its control: each group's residual permuted across genes within deciles of the reference's mean
                   expression (vec.harness.t3_decile_permutation, seed k7_wt_revert.PERM_SEED)
  resid_revert_up  the same residual on the upload's form (the central composition)
  resid_perm_up    its control on the upload's form

Worlds. Tasks g1 and g2 (WT E8.75 to Gata4 replicate 1 and 2), two replicates of one knockout run, one fold of evidence.
  plain  K8's: vec.harness.run_task on k7_wt_revert's tasks, the harness's scorer, the reference as released.
  dm     U3d's, through orth_diag.py's own functions: the frozen script's context with the two replicates
         (t3_decompose.build_context, the gata4 preset), per (replicate, seed) its truth half, reference working copy,
         per-seed estimate and per-label depth ratios; the reference working copy and the wild-type pool thinned to the
         replicate's depth (orth_diag.matched_counts on its streams); anchors recomputed on the matched reference. The
         arms draw from the matched pool. A check arm, `dm_B`, is U3d's own `dm:transfer:mab21l2_all` (the barred
         factor on the matched pool at the reference's composition, on U3d's stream) and is compared with U3d's
         committed rows: it can only block the question, never void the run. The world sets no depth in any upload;
         no T3 depth is derived for a file.

The reading (fixed in the row before any run; `decide` and `guard`):
  GUARD     the plain world's refused_form against K8b's committed scores (results/k8_s2_form/tasks): on each task the
            ten-seed mean within GUARD_MEAN (0.1) board points and every seed within GUARD_SEED (0.5); K17's tolerance.
            Else the run is void.
  QUESTION  barred_form beats refused_form by a paired mean above 0.3 board points with at least 8 of 10 seeds ahead
            on g1 AND on g2 in BOTH worlds, and barred_form's severity is scored (not gated) in at least 9 of 10 seeds
            on both tasks in the depth-matched world, and the check arm dm_B is within the guard's tolerance of U3d's
            committed rows on both replicates. Then a question to the captain to swap the last upload to the barred s1
            is put. Anything else: the refused s1 stands and no question is put.
  VOID      also unless every file holds the ten seeds, every arm, a clean tree, one commit and one machine.
  INCOMPLETE  a missing task file (a run that crashed, was refused or killed) is not a reading: `read` opens nothing,
            writes nothing and exits 3; the world is rerun at the row's commit before any task file is read.

  run   --world plain --task g1|g2 --out DIR     the plain world, one task per process
  run   --world dm --out DIR [--workers N]       the depth-matched world, both replicates in one process
  read  --runs DIR --out DIR                     apply the reading; write reading.json and reading.txt
Each run writes its task file and, beside it, <file>.prov.json: the commit, the clean flag, the machine and the
file's sha256.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import socket
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (the wild-type-only factors, the tasks, the permutation seed, _paired)
import k8_s2_form as k8  # noqa: E402  (the upload's form: the central composition and its stream)
import k16_revert_dose as k16  # noqa: E402  (move_by_group: the per-group step of K8b with the exponent a parameter)
import k17_revert_dose_reread as k17  # noqa: E402  (the guard's tolerance)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

STREAM = k8.STREAM  # "t3_s2_form": every arm of a world draws from it
REFUSED, BARRED = "refused_form", "barred_form"
RESID_K7, RESID_K7_PERM = "resid_revert_k7", "resid_perm_k7"
RESID_UP, RESID_UP_PERM = "resid_revert_up", "resid_perm_up"
RULE_ARMS = (REFUSED, BARRED)
DIAG_ARMS = (RESID_K7, RESID_K7_PERM, RESID_UP, RESID_UP_PERM)
ARMS = (*RULE_ARMS, *DIAG_ARMS)
# arm -> (form, step): form "upload" is the central composition, "k7" the reference's own; step names the factor
ARM_SPEC = {REFUSED: ("upload", "revert"), BARRED: ("upload", "barred"), RESID_K7: ("k7", "resid"),
            RESID_K7_PERM: ("k7", "resid_perm"), RESID_UP: ("upload", "resid"), RESID_UP_PERM: ("upload", "resid_perm")}
CONTROL_OF = {RESID_K7_PERM: RESID_K7, RESID_UP_PERM: RESID_UP}
WORLDS = ("plain", "dm")
WORLD_NAMES = {"plain": "K8's: the harness, the reference as released",
               "dm": "U3d's: the frozen script, the reference and the pool at each replicate's per-label depth"}
TASKS = dict(k8.TASKS)  # g1 and g2
KEYS = frozenset(t.key for t in TASKS.values())
DM_TARGETS = {"g1": "rep1", "g2": "rep2"}  # the frozen script's target names (orth_diag.context_paths)
DM_CHECK = "dm_B"  # U3d's own arm, rebuilt for the check against its committed rows
MIN_GAIN, MIN_AHEAD = k7.MIN_GAIN, k7.MIN_AHEAD  # 0.3 and 8 of 10
MIN_SCORED = 9
GUARD_MEAN, GUARD_SEED = k17.GUARD_MEAN, k17.GUARD_SEED  # 0.1 and 0.5 board points
BARRED_SLOT = ("barred", "s1")
DECISION = HERE / "results" / "gata4" / "gata4_decisions.json"
RUN_DIR = REPO / "runs" / "t3_decomp" / "gata4_20261028T122425Z"
K8_TASKS = HERE / "results" / "k8_s2_form" / "tasks"  # the guard's committed scores (K8b, k8.GROUPED)
U3D_RECORDS = HERE / "results" / "orth_diag" / "orth_records.csv"  # the dm check's committed rows
_BARRED: dict = {}  # the barred factor in the file's gene order, loaded by `run` before any candidate builds
_FOLD: dict = {}


# ---------------------------------------------------------------------------------------------------------------------
# The factors


def residual_factors(per: dict, pooled: np.ndarray, clip: tuple[float, float] = k7.CLIP) -> dict:
    """Each group's factor with its component along the pooled factor removed, in log space over the genes: r_g - b_g p
    with r_g = log(per[g]), p = log(pooled) and b_g = (r_g . p) / (p . p) (no centring; a zero p leaves r_g), returned
    as exp of that, clipped to k7's clip. A group's factor that is a power of the pooled one gives 1 everywhere."""
    p = np.log(np.asarray(pooled, dtype=np.float64))
    pp = float(p @ p)
    out = {}
    for g, v in per.items():
        r = np.log(np.asarray(v, dtype=np.float64))
        b = float(r @ p) / pp if pp > 0 else 0.0
        out[g] = np.clip(np.exp(r - b * p), *clip)
    return out


def permuted(per: dict, basis: np.ndarray) -> dict:
    """K7c's permutation applied to each group's factor: within deciles of `basis`, seed k7.PERM_SEED, one gene
    permutation for every group (as k16.permuted_factors)."""
    return {g: harness.t3_decile_permutation(v, basis, k7.PERM_SEED) for g, v in per.items()}


def barred_spec(decision: Path = DECISION) -> tuple[dict, dict, dict]:
    """(the barred s1's spec, the decision file, the rules) as the builder resolves them: the hashes checked, the spec
    from resolve_spec and validate_spec. Refuses a spec with an extra (depth, stage, coverage) or whose composition is
    not the upload's (k8.spec), since the try reads the two files on one composition."""
    bs = harness._t3_head()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    hashes = bs.check_hashes(dec, allow_mismatch=False)
    spec, _ = bs.resolve_spec(dec, rules, *BARRED_SLOT, hashes)
    bs.validate_spec(spec, rules)
    if any(spec.get(k) for k in bs.EXTRAS):
        raise RuntimeError(f"the barred s1 carries an extra step: {[k for k in bs.EXTRAS if spec.get(k)]}")
    tried, comp = k8.spec()["composition"], spec["composition"]
    if any(comp.get(k) != tried[k] for k in ("kind", "depleted_states", "log2_magnitude")):
        raise RuntimeError("the barred s1's composition is not the upload's central composition")
    return spec, dec, rules


def barred_factor(run_dir: Path, dec: dict, spec: dict, genes: list[str]) -> tuple[np.ndarray, dict]:
    """The barred s1's factor file through the builder's load_factor (sha256 against the decision file; genes in
    `genes`' order, which must be the file's)."""
    bs = harness._t3_head()
    return bs.load_factor(Path(run_dir), bs.factor_name(spec["factor"]["source"], dec["preset"]), genes, dec)


def in_gene_order(phi: np.ndarray, have: list[str], want: list[str]) -> np.ndarray:
    """phi over the genes `have`, put in the order `want` (every gene of `want` must be in `have`)."""
    pos = {g: i for i, g in enumerate(have)}
    missing = [g for g in want if g not in pos]
    if missing:
        raise RuntimeError(f"{len(missing)} genes of the reference are not in the barred factor, e.g. {missing[:5]}")
    return np.asarray(phi, dtype=np.float64)[[pos[g] for g in want]]


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """What build_arm needs for one fold: the upload's composition, the wild-type-only factors, their residuals and
    the residuals' permutation, the barred factor in the reference's gene order and the barred spec."""
    f, info = bs.composition(ref.V, k8.spec(), {})
    if info["states_missing"]:
        raise RuntimeError(f"states missing from the reference's vocabulary: {info['states_missing']}")
    res = residual_factors(per, phi)
    return {"f_upload": f, "per": per, "phi": phi, "group": group, "res": res, "res_perm": permuted(res, basis),
            "phi_barred": phi_barred, "barred_spec": spec, "ones": np.ones_like(np.asarray(phi, dtype=np.float64))}


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """One arm's drawn cells (build_submission.Built), the same code on both worlds. Every arm draws from
    stable_rng(seed, STREAM); the per-group step of every arm but the barred one runs on stable_rng(seed,
    STREAM:grouped) at k7.EXPONENT (k16.move_by_group, the step of K8b), a group without its own factor taking the
    pooled factor (refused_form) or no change (the residual arms, whose pooled residual is zero)."""
    form, step = ARM_SPEC[arm]
    rng = bs.t3d.stable_rng(seed, STREAM)
    if step == "barred":
        kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
        if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
            raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
        return bs.construct(ref, n, rng, **kw)
    built = bs.construct(ref, n, rng, f=F["f_upload"] if form == "upload" else ref.V.p)
    per, fallback = {"revert": (F["per"], F["phi"]), "resid": (F["res"], F["ones"]),
                     "resid_perm": (F["res_perm"], F["ones"])}[step]
    return k16.move_by_group(bs, ref, built, F["group"], per, fallback, k7.EXPONENT,
                             bs.t3d.stable_rng(seed, f"{STREAM}:grouped"))


# ---------------------------------------------------------------------------------------------------------------------
# The plain world (K8's: the harness)


def _plain_fold(ctx) -> dict:
    """The fold's factors, made once per prepared fold (none depends on the seed)."""
    prep = k7._prepared(ctx)
    if _FOLD.get("prep") is not prep:
        if not _BARRED:
            raise RuntimeError("the barred factor is not loaded (k18_built_files.py run loads it before the harness)")
        bs = harness._t3_head()
        ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
        F = fold_factors(bs, ref, prep["per"], prep["phi"], prep["group"], harness.t3_mean_expression(pool),
                         in_gene_order(_BARRED["phi"], _BARRED["genes"], list(ref.genes)), _BARRED["spec"])
        _FOLD.clear()
        _FOLD.update(prep=prep, F=F)
    return _FOLD["F"]


def _build_plain(ctx, arm: str):
    bs = harness._t3_head()
    prep = k7._prepared(ctx)
    F = _plain_fold(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    built = build_arm(bs, ref, arm, harness.count_for(ctx.board, pool.n_obs), ctx.seed, F)
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    def one(arm):
        form, step = ARM_SPEC[arm]
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build_plain(ctx, arm),
                                 params={"form": form, "step": step, "stream": STREAM,
                                         "composition": k8.spec()["composition"] if form == "upload" else "own",
                                         "factor": ("the barred s1's (build_submission.py --branch barred "
                                                    "--submission s1)" if step == "barred" else k7.WT_REVERT_RULE),
                                         "perm_seed": k7.PERM_SEED if step == "resid_perm" else None},
                                 requires=lambda t: t.key in KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in ARMS)]


def load_barred_for_plain(data: Path, run_dir: Path, decision: Path) -> dict:
    """The barred spec and factor, in the gene order of WT E8.75 as the builder loads it (load_reference: load_stage,
    duplicates dropped), held for the plain world's builds."""
    from vec.prep import load_stage

    spec, dec, _ = barred_spec(decision)
    genes = list(load_stage(Path(data) / "T2" / "E8.75.h5ad", drop_duplicates=True).var_names)
    phi, rec = barred_factor(run_dir, dec, spec, genes)
    _BARRED.clear()
    _BARRED.update(phi=phi, genes=genes, spec=spec, record=rec)
    return rec


# ---------------------------------------------------------------------------------------------------------------------
# The depth-matched world (U3d's: orth_diag.py's functions)


_STATE: dict = {}


def _scored(X, truth, ref, anchors: dict, up, dn, seed: int) -> dict:
    """The frozen script's scoring of one prediction against one reference, as orth_diag.score_row's `record` makes
    it, with the severity state of vec.evaluate.t3_response_detail."""
    from vec.evaluate import Stage, aggregate, panel_metrics, skills_from_anchors, t3_response_detail

    vals = panel_metrics("T3", Stage(X=X), truth, ref, seed=seed)
    skills = skills_from_anchors("T3", vals, anchors)
    score, _ = aggregate("T3", skills)
    sev = t3_response_detail(X, truth.X, ref.X, (up, dn), vals)["severity"]
    return {"score": float(score), "skills": {m: float(s) for m, s in skills.items()},
            "severity": {"state": sev["state"], "beta": sev["beta"], "r2": sev["r2"]}}


def dm_row(st: dict, ti: int, seed: int) -> dict:
    """One (replicate, seed) of the depth-matched world: orth_diag.score_row's truth half, reference working copy,
    n_pred, per-seed estimate and matched reference and pool (its streams), anchors on the matched reference; then
    U3d's check arm and every arm of this try built from the matched pool."""
    import orth_diag as od
    from common.core_metrics import de_genes

    from vec.boards import RANKED_METRICS
    from vec.evaluate import Stage, panel_metrics, split_half

    bs = harness._t3_head()
    t3d = bs.t3d
    ctx = st["ctx"]
    T, pr, V = ctx.targets[ti], ctx.rules["protocol"], ctx.V
    t0 = time.time()
    rows_t = t3d.working_copy_rows(len(T.X_rel), pr["working_copy_fraction"], seed)
    rows_r = t3d.working_copy_rows(len(ctx.ref_X), pr["working_copy_fraction"], seed + 1)
    truth, other = split_half(Stage(X=T.X_rel[rows_t]), seed=pr["split_seed"])
    ref = Stage(X=ctx.ref_X[rows_r])
    n_pred = pr["prediction_cells_multiplier"] * max(truth.n, ref.n)
    pool_rows = np.flatnonzero(~np.isin(T.hash_pool, T.hash_rel[rows_t]))
    est = t3d.estimate(ctx, T, pool_rows, t3d.stable_rng(seed, f"estimate:{T.name}"))
    wt = ctx.wt.pool
    C_ref = od.matched_counts(st["ref_counts"][rows_r], st["ref_lab"][rows_r], est.ratio, st["ref_props"][rows_r],
                              t3d.stable_rng(seed, f"dm_ref:{T.name}"))
    C_pool = od.matched_counts(wt.C, wt.lab, est.ratio, wt.P, t3d.stable_rng(seed, f"dm_pool:{T.name}"))
    ref_m = Stage(X=t3d.lognorm(C_ref))
    src = t3d.Pool(C_pool, wt.lab, wt.P, V.L)
    floor_v = panel_metrics("T3", ref_m, truth, ref_m, seed=seed)
    ceil_v = panel_metrics("T3", other, truth, ref_m, seed=seed)
    anchors = {m: {"floor": floor_v.get(m), "ceiling": ceil_v.get(m)} for m in RANKED_METRICS["T3"]
               if floor_v.get(m) is not None and ceil_v.get(m) is not None}
    up, dn, _ = de_genes(truth.X, ref_m.X)
    check = _scored(t3d.build(ctx, n_pred, t3d.stable_rng(seed, od.B_STREAM), f=V.p, src=src, phi=st["phi_b"]),
                    truth, ref_m, anchors, up, dn, seed)
    refm = bs.Reference(preset=ctx.preset_name, rules=ctx.rules, genes=list(ctx.genes), V=V, pool=src,
                        coords=np.zeros((wt.C.shape[0], 3)))
    out = {}
    for arm in ARMS:
        out[arm] = _scored(build_arm(bs, refm, arm, n_pred, seed, st["F"]).X, truth, ref_m, anchors, up, dn, seed)
    print(f"[{T.name} seed {seed}] K18 depth-matched: {len(out)} arms and the check, {time.time() - t0:.0f}s",
          flush=True)
    return {"target": T.name, "seed": int(seed), "scores": {a: v["score"] for a, v in out.items()},
            "skills": {a: v["skills"] for a, v in out.items()}, "t3": {a: {"severity": v["severity"]}
                                                                      for a, v in out.items()},
            "check": {DM_CHECK: check["score"]},
            "depth_kept": {"reference": float(C_ref.sum() / max(st["ref_counts"][rows_r].sum(), 1)),
                           "pool": float(C_pool.sum() / max(wt.C.sum(), 1))}}


def _dm_task(task) -> dict:
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=_STATE["threads"]):
        return dm_row(_STATE, *task)


def run_dm(data: Path, run_dir: Path, decision: Path, workers: int, threads: int) -> dict:
    """The depth-matched world over both replicates and the ten seeds, set up as orth_diag.cmd_score sets up U3d:
    the frozen script's context with the replicates, its self-check, the barred factor through load_factor, the
    reference's counts; plus the wild-type-only factors through the builder's wt_revert_factors."""
    import multiprocessing as mp

    import orth_diag as od

    bs = harness._t3_head()
    t3d = bs.t3d
    spec, dec, rules = barred_spec(decision)
    if list(rules["protocol"]["seeds"]) != list(SEEDS):
        raise RuntimeError("the frozen script's seeds are not vec.adoption.SEEDS")
    paths = od.context_paths(Path(data), with_targets=True)
    n_par = workers * threads
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=threads, out=Path(tmp))
        if sorted(T.name for T in ctx.targets) != sorted(DM_TARGETS.values()):
            raise RuntimeError(f"targets {[T.name for T in ctx.targets]}, expected {sorted(DM_TARGETS.values())}")
        t3d.selfcheck(ctx, SEEDS[0])
        genes = list(ctx.genes)
        phi_b, frec = barred_factor(run_dir, dec, spec, genes)
        ref_C, ref_lab, ref_P, gap = od.reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        wt_path = Path(paths["wt"])
        series = bs.wt_revert_series(wt_path, {n: Path(data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                               if n != wt_path.stem})
        phi, per, group, rrec = bs.wt_revert_factors(wt_path, series, genes, ctx.wt.pool.C)
        refw = bs.Reference(preset=ctx.preset_name, rules=rules, genes=genes, V=ctx.V, pool=ctx.wt.pool,
                            coords=np.zeros((ctx.wt.pool.C.shape[0], 3)))
        F = fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b, spec)
        _STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F, threads=threads)
        tasks = [(ti, s) for ti in range(len(ctx.targets)) for s in SEEDS]
        rows: list[dict] = []
        if workers > 1:
            if not bs.openmp_limited():
                raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
            with mp.get_context("fork").Pool(min(workers, len(tasks))) as p:
                rows.extend(p.imap_unordered(_dm_task, tasks))
        else:
            rows.extend(_dm_task(t) for t in tasks)
    by_name = {v: k for k, v in DM_TARGETS.items()}
    files = [*(Path(data) / f for f in od.WT_FILES.values()), *(Path(data) / f for f in od.KO_FILES.values()),
             *series.values(), Path(decision), Path(run_dir) / "factors" / f"{od.BARRED}.csv"]
    return {"world": "dm", "seeds": list(SEEDS), "arms": list(ARMS),
            "rows": {by_name[n]: sorted((r for r in rows if r["target"] == n), key=lambda r: r["seed"])
                     for n in DM_TARGETS.values()},
            "barred_factor": frec, "wt_revert": {k: v for k, v in rrec.items() if k != "series"},
            "reference_counts_largest_gap": gap,
            "files_sha256": {str(f): t3d.sha256(f) for f in dict.fromkeys(files)}}


# ---------------------------------------------------------------------------------------------------------------------
# Running


def stamp() -> dict:
    """The commit, the clean flag (vec.harness's rule: tracked files only), the machine and the time."""
    return {"commit": harness._git("rev-parse", "HEAD"),
            "tree_clean": harness._git("status", "--porcelain", "--untracked-files=no") == "",
            "machine": socket.gethostname(), "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}


def write_task_file(path: Path, res: dict) -> dict:
    """The task file and, beside it, <file>.prov.json with the commit, the clean flag, the machine and its sha256.
    The task file is written to a temporary name and renamed, so a run killed while writing leaves no task file and
    the reading sees it as missing (not a reading), never as a short file."""
    raw = json.dumps(res, indent=1, default=float).encode()
    tmp = path.parent / f".{path.name}.part"
    tmp.write_bytes(raw)
    tmp.replace(path)
    prov = {k: res["provenance"][k] for k in ("commit", "tree_clean", "machine")} | {
        "file": path.name, "sha256": hashlib.sha256(raw).hexdigest()}
    (path.parent / f"{path.name}.prov.json").write_text(json.dumps(prov, indent=1) + "\n")
    return prov


def cmd_run(a) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    name = f"plain_{a.task}.json" if a.world == "plain" else "dm.json"
    path = out / name
    if path.exists():
        sys.exit(f"{path} exists")
    if a.world == "plain":
        if a.task not in TASKS:
            sys.exit("--world plain needs --task g1 or g2")
        st = stamp()
        rec = load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
        res = harness.run_task(TASKS[a.task], candidates(), seeds=SEEDS, data_dir=Path(a.data))
        res["provenance"]["machine"] = st["machine"]
        res |= {"world": "plain", "task_source": "experiments/t3_decomp/k18_built_files.py",
                "factor_record": k7._PREP["fold"]["record"], "barred_factor": rec}
    else:
        st = stamp()
        res = run_dm(Path(a.data), Path(a.run_dir), Path(a.decision), a.workers, a.threads)
        res |= {"task_source": "experiments/t3_decomp/k18_built_files.py", "provenance": st}
    prov = write_task_file(path, res)
    print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}, sha256 "
          f"{prov['sha256']}", flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def distance(rows: list[dict], arm: str, ref_rows: list[dict], ref_arm: str) -> dict | None:
    """`arm` against the committed `ref_arm`, seed by seed: the difference of the ten-seed means and the largest
    absolute per-seed difference (None when a seed is missing from either)."""
    ref = {x["seed"]: x["scores"].get(ref_arm) for x in ref_rows}
    if not rows or any(ref.get(q["seed"]) is None or arm not in q["scores"] for q in rows):
        return None
    d = np.array([q["scores"][arm] - ref[q["seed"]] for q in rows], dtype=float)
    return {"mean_diff": float(d.mean()), "max_abs_diff": float(np.abs(d).max()),
            "seeds_off": [q["seed"] for q, v in zip(rows, d, strict=True) if abs(v) > 1e-9]}


def guard_passes(d: dict | None) -> bool:
    """K17's tolerance: the means within GUARD_MEAN and every seed within GUARD_SEED."""
    return d is not None and abs(d["mean_diff"]) <= GUARD_MEAN and d["max_abs_diff"] <= GUARD_SEED


def committed_k8(name: str, refs: Path = K8_TASKS) -> list[dict]:
    return json.loads((refs / f"{name}.json").read_text())["rows"]


def _median(v: list) -> float | None:
    v = [float(x) for x in v if x is not None]
    return float(np.median(v)) if v else None


def summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """Per arm: the board score's mean and sd, each metric's mean skill, the severity states, the seeds scored and
    the median beta and r2; the paired differences against refused_form (and, where present, a control against its
    arm and refused_form against the floor)."""
    out: dict = {"arms": {}, "paired": {}}
    for x in arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "severity_states": {st: sum(v["state"] == st for v in sev) for st in sorted({v["state"] for v in sev})},
            "severity_scored": int(sum(v["state"] == "scored" for v in sev)),
            "median_beta": _median([v.get("beta") for v in sev]), "median_r2": _median([v.get("r2") for v in sev])}
    pairs = [(x, REFUSED) for x in arms if x not in (REFUSED, harness.FLOOR)]
    pairs += [(c, a) for c, a in CONTROL_OF.items() if c in arms and a in arms]
    if harness.FLOOR in arms:
        pairs.append((REFUSED, harness.FLOOR))
    out["paired"] = {f"{x} - {y}": k7._paired(rows, x, y) for x, y in pairs}
    return out


def decide(per: dict, check: dict | None) -> dict:
    """The rule fixed in row K18, on summarise()'s output per world and task and on dm_check()'s output. The check
    arm can only block the question: outside the guard's tolerance on either replicate (or not computable), the
    depth-matched world's setup is not shown to reproduce U3d, and no question is put; the run is not void.
    refused_form's scored count in the depth-matched world is reported beside barred_form's and decides nothing."""
    key = f"{BARRED} - {REFUSED}"
    beats = {w: {n: bool(per[w][n]["paired"][key]["mean"] > MIN_GAIN and per[w][n]["paired"][key]["ahead"] >= MIN_AHEAD)
                 for n in TASKS} for w in WORLDS}
    scored = {n: per["dm"][n]["arms"][BARRED]["severity_scored"] for n in TASKS}
    scored_refused = {n: per["dm"][n]["arms"][REFUSED]["severity_scored"] for n in TASKS}
    scored_ok = all(v >= MIN_SCORED for v in scored.values())
    check_ok = bool(check) and all(check.get(n, {}).get("within_tolerance") for n in TASKS)
    question = all(beats[w][n] for w in WORLDS for n in TASKS) and scored_ok and check_ok
    return {"beats": beats, "dm_barred_severity_scored": scored, "dm_refused_severity_scored": scored_refused,
            "dm_scored_ok": bool(scored_ok), "dm_check_ok": check_ok, "question": bool(question),
            "outcome": ("a question to the captain: swap T3's last upload to the barred s1 (a priced question after a "
                        "cold check; the try changes no entry)" if question else
                        "no question; the refused s1 stands as T3's last upload")}


def _read_file(p: Path, why: list[str]) -> tuple[dict | None, str | None]:
    if not p.exists():
        why.append(f"missing {p.name}")
        return None, None
    raw = p.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def _rows_ok(name: str, rows: list[dict], arms: tuple[str, ...], why: list[str]) -> bool:
    """The ten seeds, and every arm in every row's scores, skills and severity (the plain world's arms include the
    floor), so that a short row voids the reading instead of raising in summarise."""
    if [x["seed"] for x in rows] != list(SEEDS):
        why.append(f"{name}: other seeds")
        return False
    miss = sorted({a for q in rows for a in arms
                   if a not in q.get("scores", {}) or a not in q.get("skills", {}) or a not in q.get("t3", {})})
    if miss:
        why.append(f"{name}: missing arms {miss}")
        return False
    return True


def dm_check(rows: dict, records: Path = U3D_RECORDS) -> dict:
    """The dm world's check arm against U3d's committed `dm:transfer:mab21l2_all` rows, per replicate: distance()'s
    figures and whether they are within the guard's tolerances. Outside them on either replicate, decide() puts no
    question; the run is not void (K16's lesson: a reproduction check across machines must not void a try)."""
    import pandas as pd

    rec = pd.read_csv(records)
    rec = rec[(rec["world"] == "depth_matched") & (rec["variant"] == "dm:transfer:mab21l2_all")]
    out = {}
    for n, target in DM_TARGETS.items():
        ref_rows = [{"seed": int(r.seed), "scores": {DM_CHECK: float(r.t3)}} for r in rec[rec["target"] == target]
                    .itertuples()]
        mine = [{"seed": q["seed"], "scores": {DM_CHECK: q["check"][DM_CHECK]}} for q in rows.get(n, [])]
        d = distance(mine, DM_CHECK, ref_rows, DM_CHECK)
        out[n] = {"distance": d, "within_tolerance": guard_passes(d)}
    return out


TASK_FILES = (*(f"plain_{n}.json" for n in TASKS), "dm.json")
INCOMPLETE = ("not a reading: a run ended without writing its task file; nothing was read. Row K18: rerun that "
              "world at the row's commit (after a committed code fix noted as an amendment to the row, if the cause "
              "is in the code) before any task file is read")


def reading(runs: Path, refs: Path = K8_TASKS, u3d: Path = U3D_RECORDS) -> dict:
    """Apply the reading of row K18 to a run folder (plain_g1.json, plain_g2.json, dm.json). A run that crashed, was
    refused or was killed writes no task file; while any of the three is missing, nothing in the folder is opened and
    the result is 'incomplete', not void, so the missing world can be rerun blind."""
    missing = [f for f in TASK_FILES if not (runs / f).exists()]
    if missing:
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    why, files, per = [], {}, {"plain": {}, "dm": {}}
    prov: dict = {}
    for n in TASKS:
        r, h = _read_file(runs / f"plain_{n}.json", why)
        if r is None:
            continue
        files[f"plain_{n}.json"] = h
        prov[f"plain_{n}.json"] = {k: r["provenance"].get(k) for k in ("commit", "tree_clean", "machine")}
        rows = sorted(r["rows"], key=lambda x: x["seed"])
        if r["task"]["key"] != TASKS[n].key:
            why.append(f"plain_{n}: another task")
            continue
        if _rows_ok(f"plain_{n}", rows, (harness.FLOOR, *ARMS), why):
            per["plain"][n] = summarise(rows, (harness.FLOOR, *ARMS)) | {
                "guard": distance(rows, REFUSED, committed_k8(n, refs), k8.GROUPED)}
    r, h = _read_file(runs / "dm.json", why)
    check = None
    if r is not None:
        files["dm.json"] = h
        prov["dm.json"] = {k: r["provenance"].get(k) for k in ("commit", "tree_clean", "machine")}
        for n in TASKS:
            rows = sorted(r["rows"].get(n, []), key=lambda x: x["seed"])
            if not rows:
                why.append(f"dm: no rows for {n}")
            elif _rows_ok(f"dm {n}", rows, ARMS, why):
                per["dm"][n] = summarise(rows, ARMS)
        check = dm_check(r["rows"], u3d)
    for f, p in prov.items():
        if not p["tree_clean"]:
            why.append(f"{f}: a dirty tree")
    commits = {p["commit"] for p in prov.values()}
    machines = {p["machine"] for p in prov.values()}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(map(str, commits))}")
    if len(machines) > 1:
        why.append(f"results from different machines {sorted(map(str, machines))}")
    for n, v in per["plain"].items():
        if not guard_passes(v["guard"]):
            d = v["guard"]
            why.append(f"guard: plain {n}: {REFUSED} is off K8b's committed scores" + (
                f" by {d['mean_diff']:+.3f} in the mean and {d['max_abs_diff']:.3f} in one seed" if d else
                " (a seed is missing)") + f" (limits {GUARD_MEAN} and {GUARD_SEED})")
    base = {"incomplete": False, "task_sha256": files, "provenance": prov, "dm_check": check,
            "guard": {n: v["guard"] for n, v in per["plain"].items()}}
    if why:
        return base | {"void": True, "why": why, "worlds": per,
                       "outcome": "void: no question; the refused s1 stands as T3's last upload"}
    return base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                   "worlds": per} | decide(per, check)


def _f(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def _dist(d: dict | None) -> str:
    if d is None:
        return "n/a"
    return (f"mean {d['mean_diff']:+.4f}, largest seed {d['max_abs_diff']:.4f}, {len(d['seeds_off'])} of "
            f"{len(SEEDS)} seeds differ")


def _world_lines(world: str, per: dict) -> list[str]:
    lines = []
    for n, t in TASKS.items():
        v = per.get(n)
        if v is None:
            lines.append(f"{world} {n} {t.key}: not read")
            continue
        lines.append(f"{world} {n} {t.key}:")
        for a, s in v["arms"].items():
            lines.append(f"    {a:<16} {s['mean']:.2f} (sd {s['sd']:.2f}); severity scored {s['severity_scored']} of "
                         f"{len(SEEDS)} {s['severity_states']}, median beta {_f(s['median_beta'])}, r2 "
                         f"{_f(s['median_r2'])}; skills " + ", ".join(f"{m} {x:.3f}" for m, x in s["skills"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
    return lines


def text(res: dict) -> str:
    lines = [f"Try K18: {res['outcome']}."]
    if res.get("incomplete"):
        return "\n".join([*lines, "missing task files: " + ", ".join(res["missing"])]) + "\n"
    if res["void"]:
        lines.append("VOID: " + "; ".join(res["why"]))
    else:
        lines.append(f"commit {res['commit']}, machine {res['machine']}. Barred over refused, paired mean > {MIN_GAIN} "
                     f"with >= {MIN_AHEAD} of {len(SEEDS)} ahead: " + "; ".join(
                         f"{w} {n} {'yes' if res['beats'][w][n] else 'no'}" for w in WORLDS for n in TASKS)
                     + "; barred severity scored in the depth-matched world: " + ", ".join(
                         f"{n} {c} of {len(SEEDS)}" for n, c in res["dm_barred_severity_scored"].items())
                     + f" (needs >= {MIN_SCORED} on both; refused_form, reported only: " + ", ".join(
                         f"{n} {c} of {len(SEEDS)}" for n, c in res["dm_refused_severity_scored"].items())
                     + f"); the depth-matched check within tolerance on both replicates: "
                     f"{'yes' if res['dm_check_ok'] else 'no'} (needed for a question).")
    lines.append("guard (plain refused_form against K8b's committed scores): " + "; ".join(
        f"{n} {_dist(d)}" for n, d in res["guard"].items()))
    if res.get("dm_check"):
        lines.append("check (the dm world's U3d arm against U3d's committed rows; it can only block the question): "
                     + "; ".join(
            f"{n} {_dist(c['distance'])}, within the guard's tolerance {'yes' if c['within_tolerance'] else 'no'}"
            for n, c in res["dm_check"].items()))
    lines.append("task files read: " + "; ".join(
        f"{f} sha256 {h}, commit {res['provenance'][f]['commit']}, clean {res['provenance'][f]['tree_clean']}, "
        f"machine {res['provenance'][f]['machine']}" for f, h in res["task_sha256"].items()))
    lines.append("Diagnostic arms (resid_*) decide nothing; g1 and g2 are two replicates of one knockout run.")
    for w in WORLDS:
        lines += ["", f"WORLD {w} ({WORLD_NAMES[w]}):", *_world_lines(w, res["worlds"].get(w, {}))]
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    res = reading(Path(a.runs))
    if res.get("incomplete"):  # nothing was read; nothing is written, so the rerun stays blind
        print(text(res), end="")
        sys.exit(3)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--world", choices=WORLDS, required=True)
    r.add_argument("--task", choices=sorted(TASKS))
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--run-dir", default=str(RUN_DIR))
    r.add_argument("--decision", default=str(DECISION))
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--threads", type=int, default=8)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
