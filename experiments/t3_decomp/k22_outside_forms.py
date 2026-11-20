#!/usr/bin/env python
"""Try K22 (experiments/harness/TRIES.md; the row is placed before any run, and `pregate` and `run` refuse without it):
the barred factor times one outside-knockout table, read against the barred s1 on the two Gata4 replicates in K18's two
practice worlds. A labelled reading that may name at most two T3 candidate files by the rule fixed in the row; not the
adoption rule.

Built on k18_built_files.py and in the manner of k21_factor_forms.py, both imported unchanged: K18's stream
(`t3_s2_form`), tasks (g1, g2), worlds (plain, dm), seeds, barred spec and factor, plain build, depth-matched row,
stamp, task-file writer, folder rule, masking and guard; K21's seed parser, summary and K18-row reader. For a run, three
of k18's module names are replaced inside this process only (`k18_arms`): ARMS, build_arm (whose barred_form branch is
K18's own) and fold_factors (K18's, plus the outside tables in the world's gene order and their permutations). So every
arm draws the barred s1's cells exactly as K18's barred_form did, and K18's barred_form is this try's base and guard.

The outside tables come only through experiments/t3_decomp/outside_forms.py, the module the builder
build_outside_form.py applies: by name, from the folder of the environment variable VEC_T3_EXT_TABLES, at the sha256
of the committed experiments/t3_ext/TABLES.sha256 and the pins in outside_forms.py. No path is taken from the command
line.

Arms (phi_M = transfer_mab21l2_all as loaded; T_all, T_unrel the pooled tables `all` and `unrelated`; perm =
outside_forms.permuted within deciles of the basis K18's fold_factors receives in each world, seed 20261113, 10 bins):
  barred_form    phi_M (base and guard; K18's arm value for value)
  x_all          phi_M * T_all            x_all_perm    phi_M * perm(T_all)     (control)
  x_unrel        phi_M * T_unrel          x_unrel_perm  phi_M * perm(T_unrel)   (control)
  rep_all        T_all alone, in place of phi_M (the labelled control of the planner's doubt, decision log 2026-11-20
                 (m); counted, never named)
The pre-gate (outside_forms.pregate; target-free) drops a try with its control (and rep_all with x_all); a dropped arm
is not run and not counted. When both tries drop, `pregate` writes pregate.json and exits 5: nothing runs or counts.

  pregate --out DIR                                   the pre-gate: DIR/pregate.json and its lines (no score)
  run --world plain --task g1|g2 [--seeds 0-1] --out DIR   one plain task, seed by seed (resumable)
  run --world dm [--seeds ...] [--workers 8 --threads 2] --out DIR   both replicates, row by row (resumable)
  read --runs DIR --out DIR                           the fixed reading: reading.json and reading.txt (never over an
                                                      existing one: exit 4)

Exit codes: 0 done; 1 a crash (masked traceback; seeds written so far stay); 2 refused before anything loads or when a
table loads (no row naming the script, a dirty tree, pins not filled, TABLES.sha256 untracked or off its pin, the
folder variable unset or refused, bytes or records that disagree, genes off the barred factor's, the frozen T3 files
off their pins); 3 INCOMPLETE (read: a file missing; nothing written); 4 read: the reading files exist already; 5
pregate: both tries dropped, nothing runs.

The reading (fixed in the row; `reading`, `decide`): as row K21, with these arms; VOID also unless every task file and
pregate.json record the same tables (listing sha256, table sha256, COMBINE.json sha256) and those equal the pins;
QUALIFIED (x_all, x_unrel) as K21 with (d) needing at least one dm cell where the try is AHEAD and beats its control;
NAMES: rank by AHEAD cells, then the mean of the four paired means (ties: x_all, x_unrel); file 1 the first on base B,
file 2 the second only when the pre-gate's corr(ln T_all, ln T_unrel) is below 0.9. The doubt statement (report only):
"the planner's doubt holds on this fold" when rep_all is BEHIND barred_form in all four cells, else "the doubt is not
shown on this fold".
"""

from __future__ import annotations

import argparse
import hashlib
import json
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
import k7_wt_revert as k7  # noqa: E402  (the permutation seed, MIN_GAIN and MIN_AHEAD)
import k8_s2_form as k8  # noqa: E402  (the upload's central composition, for the record)
import k18_built_files as k18  # noqa: E402  (the worlds, tasks, stream, seeds, build path and scoring; unchanged)
import k21_factor_forms as k21  # noqa: E402  (the seed parser, the summary and K18's committed rows; unchanged)
import outside_forms as ext  # noqa: E402  (the tables and the product, shared with build_outside_form.py)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

SCRIPT = "experiments/t3_decomp/k22_outside_forms.py"
ROW_SCRIPT = f"`{SCRIPT}`"  # the ledger row names this script
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"
STREAM = k18.STREAM  # "t3_s2_form"
BASE = k18.BARRED  # "barred_form"
DOUBT = "rep_all"
ARMS = (BASE, "x_all", "x_all_perm", "x_unrel", "x_unrel_perm", DOUBT)
NAMEABLE = ("x_all", "x_unrel")  # also the tie order
CONTROL_OF = {"x_all": "x_all_perm", "x_unrel": "x_unrel_perm"}
TABLE_OF = {"x_all": "all", "x_all_perm": "all", "x_unrel": "unrelated", "x_unrel_perm": "unrelated", DOUBT: "all"}
ARMS_OF_TABLE = {"all": ("x_all", "x_all_perm", DOUBT), "unrelated": ("x_unrel", "x_unrel_perm")}
TASKS = k18.TASKS
WORLDS = k18.WORLDS
MIN_GAIN, MIN_AHEAD = k7.MIN_GAIN, k7.MIN_AHEAD  # 0.3 and 8 of 10
MAX_EXTRA_GATED = 1
GUARD_MEAN, GUARD_SEED = k18.GUARD_MEAN, k18.GUARD_SEED  # 0.1 and 0.5
K18_TASKS = k21.K18_TASKS
K18_SHA256 = dict(k21.K18_SHA256)
RESULTS = HERE / "results" / "k22_outside_forms"
EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_NOTHING_RUNS = 2, 3, 4, 5
_EXT: dict = {}  # the loaded tables in the barred factor's file order, set before any arm builds


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


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


def frozen_ok() -> dict:
    """The three frozen T3 files at their pinned bytes (as the thin builders pin them)."""
    bs = harness._t3_head()
    got = {"t3_decompose.py": bs.sha256(bs.T3D_PATH), "rules.yaml": bs.sha256(bs.RULES_PATH),
           "build_submission.py": bs.sha256(bs.BUILDER_PATH)}
    bad = [k for k in ext.FROZEN_SHA256 if got.get(k) != ext.FROZEN_SHA256[k]]
    if bad:
        refuse(f"frozen file(s) {bad} differ from the pinned bytes")
    return got


DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing written. Row K22: commit or restore first; INCOMPLETE, not VOID")


def preconditions() -> tuple[str, dict]:
    """Before anything loads: the row, a clean tree, the pins filled, the listing tracked and at its pin, the frozen
    files. Returns (row id, stamp)."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    try:
        ext.read_listing()
    except ext.TableRefused as e:
        refuse(str(e))
    frozen_ok()
    return rid, st


# ---------------------------------------------------------------------------------------------------------------------
# The tables and the arms


def load_ext(genes: list[str], phi_m: np.ndarray, basis=None, den=None) -> dict:
    """The listed tables against the barred factor's file genes (order included), the pre-gate and the arms that run;
    held in _EXT for fold_factors. A refusal of outside_forms exits 2."""
    try:
        loaded = ext.load_listed(list(genes))
        pre = ext.pregate(phi_m, loaded, basis, den)
    except (ext.TableRefused, ext.FormRefused) as e:
        refuse(str(e))
    _EXT.clear()
    _EXT.update(genes=list(genes), T={n: T for n, (T, _) in loaded["tables"].items()}, loaded=loaded, pre=pre,
                record=ext.ext_record(loaded), table_records={n: r for n, (_, r) in loaded["tables"].items()})
    return pre


def dropped_arms(pre: dict) -> dict[str, str]:
    """The pre-gate's table drops as arm drops: a dropped table drops its try, its control and (all) rep_all."""
    return {a: f"table {t}: {why}" for t, why in pre["dropped"].items() for a in ARMS_OF_TABLE[t]}


def active_arms(pre: dict) -> tuple[str, ...]:
    d = dropped_arms(pre)
    return tuple(a for a in ARMS if a not in d)


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K18's fold_factors (the original, called unchanged) plus each loaded table in the reference's gene order (the
    in_gene_order K18 applies to phi_M) and its permutation within deciles of the `basis` this function receives."""
    if not _EXT:
        raise RuntimeError("the outside tables are not loaded (k22_outside_forms.py loads them before any arm)")
    F = _K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    F["ext"] = {n: k18.in_gene_order(T, _EXT["genes"], list(ref.genes)) for n, T in _EXT["T"].items()}
    F["ext_perm"] = {n: ext.permuted(T, basis) for n, T in F["ext"].items()}
    return F


def arm_factor(arm: str, F: dict) -> np.ndarray:
    """The factor an arm passes in place of phi, in the reference's gene order."""
    if arm == BASE:
        return F["phi_barred"]
    if arm == DOUBT:
        return F["ext"]["all"]
    src = F["ext_perm"] if arm in CONTROL_OF.values() else F["ext"]
    return ext.product(F["phi_barred"], 1.0, src[TABLE_OF[arm]])


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """One arm's drawn cells. Every arm draws from stable_rng(seed, STREAM) through bs.construct with the barred
    spec's keywords (prediction_inputs: the central composition, w = 1, no depth ratio or stage mixing), so every arm
    holds the barred s1's rows; only the factor passed differs. barred_form passes the keywords unchanged (K18's
    arm)."""
    if arm not in ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    if arm != BASE:
        kw["phi"] = arm_factor(arm, F)
    return bs.construct(ref, n, rng, **kw)


_K18_ORIGINAL = {"ARMS": k18.ARMS, "build_arm": k18.build_arm, "fold_factors": k18.fold_factors}


@contextmanager
def k18_arms(arms: tuple[str, ...]):
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's, restored on exit."""
    if not set(arms) <= set(ARMS) or BASE not in arms:
        raise RuntimeError(f"arms {arms} are not a set of this try's arms with the base")
    k18.ARMS, k18.build_arm, k18.fold_factors = tuple(arms), build_arm, fold_factors
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
                                                    f"outside_forms ({arm}; table {TABLE_OF[arm]})"),
                                         "perm_seed": ext.PERM_SEED if arm in CONTROL_OF.values() else None},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


# ---------------------------------------------------------------------------------------------------------------------
# Running


def progress_log(world: str, task: str):
    """The harness's per-seed line lists every arm's score; this prints the world, the task and the seed only."""
    import re

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"K22 {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


class _Unclean(RuntimeError):
    pass


def _stamp_ok(st: dict) -> None:
    now = k18.stamp()
    if now["tree_clean"] is not True or now["commit"] != st["commit"]:
        raise _Unclean(f"the tree is now commit {now['commit']}, clean {now['tree_clean']}, against {st['commit']} "
                       "clean at the start: no task file written (INCOMPLETE; rerun on a clean tree at the commit)")


def _pregate_record(pre: dict, arms: tuple[str, ...]) -> dict:
    return {"dropped": dropped_arms(pre), "arms": list(arms), "ext": dict(_EXT["record"])}


def _world_figures(F: dict) -> dict:
    """Per world, target-free: genes changed against phi_M by each try and each control in this world's order."""
    out = {}
    for n, T in F["ext"].items():
        out[n] = {"try": ext.figures(F["phi_barred"], ext.product(F["phi_barred"], 1.0, T)),
                  "control": ext.figures(F["phi_barred"], ext.product(F["phi_barred"], 1.0, F["ext_perm"][n]))}
    return out


def run_plain(a, st: dict, out: Path) -> None:
    folder = out / f"plain_{a.task}"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    pre = load_ext(k18._BARRED["genes"], k18._BARRED["phi"])
    arms = active_arms(pre)
    if arms == (BASE,):
        refuse("the pre-gate drops both tries: nothing runs")
    with k18_arms(arms):
        for s in k21.parse_seeds(a.seeds):
            path = folder / f"seed_{s}.json"
            if path.exists():
                print(f"K22 plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            res = harness.run_task(TASKS[a.task], candidates(arms), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log("plain", a.task))
            got = res["provenance"]
            if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
                raise _Unclean(f"the harness's stamp says commit {got.get('commit')}, clean {got.get('tree_clean')}: "
                               "no task file written (INCOMPLETE)")
            _stamp_ok(st)
            res["provenance"]["machine"] = st["machine"]
            fold = k18._FOLD.get("F")
            res |= {"world": "plain", "task_source": SCRIPT, "seed": s, "pregate": _pregate_record(pre, arms),
                    "factor_record": k7._PREP["fold"]["record"], "barred_factor": rec,
                    "world_figures": _world_figures(fold) if fold else None}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)


def run_dm(a, st: dict, out: Path) -> None:
    """The depth-matched world, set up exactly as k18.run_dm and k21.run_dm set it up, then each (replicate, seed) row
    through k18.dm_row with this try's arms; rows already written at this commit on this machine are skipped."""
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
        pre = load_ext(genes, phi_b)
        arms = active_arms(pre)
        if arms == (BASE,):
            refuse("the pre-gate drops both tries: nothing runs")
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
            wf = _world_figures(F)
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder / f"{by_name[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                _stamp_ok(st)
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


def pregate_lines(pre: dict) -> list[str]:
    lines = [f"pre-gate (the tables and phi_M alone; no score): {pre['n_genes']} genes; phi_M "
             + ", ".join(f"{k} {_fmt(v)}" for k, v in pre["phi_m"].items())]
    for kind in ("tables", "products", "controls"):
        for n, r in pre[kind].items():
            lines.append(f"  {kind[:-1]} {n}: " + ", ".join(f"{k} {_fmt(v)}" for k, v in r.items()))
    lines.append(f"  corr(ln T_all, ln T_unrel) over genes != 1 in either: {_fmt(pre['corr_ln_all_ln_unrelated'])}")
    lines.append("  dropped: " + ("; ".join(f"{k} ({v})" for k, v in pre["dropped"].items()) or "none"))
    return lines


def cmd_pregate(a) -> None:
    """The pre-gate before any scoring: phi_M through the builder's loader (WT E8.75's gene order), the tables through
    outside_forms, the controls' figures on WT E8.75's mean expression and the den-weighted totals on its mean
    within-cell proportions (figures only; no drop depends on them); DIR/pregate.json."""
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
        pre = load_ext(genes, phi_m, harness.t3_mean_expression(wt), den)
        arms = active_arms(pre)
        res = {"pregate": pre, "arms": list(arms), "dropped_arms": dropped_arms(pre), "row": rid,
               "barred_factor": rec, "ext": _EXT["record"], "table_records": _EXT["table_records"],
               "provenance": st}
        k18.write_task_file(path, res)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print("\n".join(pregate_lines(pre)), flush=True)
    if arms == (BASE,):
        print("both tries dropped: nothing runs and nothing is counted", flush=True)
        sys.exit(EXIT_NOTHING_RUNS)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _cells():
    return [(w, n) for w in WORLDS for n in TASKS]


def _summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """k21.summarise with this try's control pairs."""
    saved = dict(k21.CONTROL_OF)
    k21.CONTROL_OF.clear()
    k21.CONTROL_OF.update(CONTROL_OF)
    try:
        return k21.summarise(rows, arms)
    finally:
        k21.CONTROL_OF.clear()
        k21.CONTROL_OF.update(saved)


def _beats(v: dict) -> bool:
    return bool(v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD)


def doubt_statement(per: dict, arms: tuple[str, ...]) -> dict:
    if DOUBT not in arms:
        return {"ran": False, "statement": "rep_all was not run (dropped with x_all)"}
    means = {f"{w} {n}": per[w][n]["paired"][f"{DOUBT} - {BASE}"]["mean"] for w, n in _cells()}
    behind = [k for k, m in means.items() if m < -MIN_GAIN]
    return {"ran": True, "behind_cells": behind, "means": means,
            "statement": ("the planner's doubt holds on this fold" if len(behind) == 4 else
                          "the doubt is not shown on this fold")}


def decide(per: dict, arms: tuple[str, ...], corr_all_unrel: float | None) -> dict:
    """The rule fixed in the row, on the summaries per (world, task)."""
    status: dict = {}
    for t in NAMEABLE:
        if t not in arms:
            status[t] = {"ran": False, "qualified": False}
            continue
        c = CONTROL_OF[t]
        d = {(w, n): per[w][n]["paired"][f"{t} - {BASE}"] for w, n in _cells()}
        ahead = {k: _beats(v) for k, v in d.items()}
        behind = {k: bool(v["mean"] < -MIN_GAIN) for k, v in d.items()}
        gated = {k: (per[k[0]][k[1]]["arms"][t]["severity_gated"], per[k[0]][k[1]]["arms"][BASE]["severity_gated"])
                 for k in d}
        ctrl = {k: per[k[0]][k[1]]["paired"][f"{t} - {c}"] for k in d}
        a = (any(ahead[("plain", n)] for n in TASKS) and any(ahead[("dm", n)] for n in TASKS)
             and any(ahead[(w, "g1")] for w in WORLDS) and any(ahead[(w, "g2")] for w in WORLDS))
        b = not any(behind.values())
        cc = all(g <= g0 + MAX_EXTRA_GATED for g, g0 in gated.values())
        dd = (all(_beats(ctrl[k]) for k in d if ahead[k])
              and any(ahead[("dm", n)] and _beats(ctrl[("dm", n)]) for n in TASKS))
        status[t] = {"ran": True, "qualified": bool(a and b and cc and dd), "a": bool(a), "b": bool(b), "c": bool(cc),
                     "d": bool(dd), "ahead_cells": [f"{w} {n}" for (w, n), v in ahead.items() if v],
                     "behind_cells": [f"{w} {n}" for (w, n), v in behind.items() if v],
                     "mean_of_four": float(np.mean([v["mean"] for v in d.values()])),
                     "gated_vs_base": {f"{w} {n}": list(v) for (w, n), v in gated.items()}}
    q = [t for t in NAMEABLE if status[t]["qualified"]]
    ranked = sorted(q, key=lambda t: (-len(status[t]["ahead_cells"]), -status[t]["mean_of_four"], NAMEABLE.index(t)))
    file2 = len(ranked) > 1 and corr_all_unrel is not None and corr_all_unrel < ext.FILE2_MAX_CORR
    names = ranked[:2] if file2 else ranked[:1]
    kinds = {"x_all": "xall", "x_unrel": "xunrel"}
    return {"status": status, "ranked": ranked, "names": names,
            "file2_rule": {"corr_ln_all_ln_unrelated": corr_all_unrel, "below": ext.FILE2_MAX_CORR,
                           "applies": len(ranked) > 1, "named": bool(file2)},
            "files": [f"{kinds[t]}-B (base B of decision log 2026-11-20 (q); build_outside_form.py --form "
                      f"{kinds[t]}-<B>; table {TABLE_OF[t]})" for t in names],
            "outcome": ("names " + " and ".join(f"file {i + 1}: {t} on base B" for i, t in enumerate(names))
                        + "; each after the C list and after any K21 file, under the R14/R17 conditions" if names
                        else "no file: nothing qualifies; T3's list is unchanged")}


def expected_files() -> list[str]:
    return k21.expected_files()


INCOMPLETE = ("not a reading: pregate.json or a seed's task file is missing (a crash, a refusal such as a dirty tree, "
              "out of memory or a kill); nothing was read. Row K22: if the cause is outside the code, rerun the "
              "missing seeds into this folder at the run's commit on the same machine (the runs skip what is there); "
              "after a committed code fix, rerun everything at the new commit into a fresh folder and never read this "
              "one")


def pinned_ext() -> dict:
    return {"listing_sha256": ext.K22_LISTING_SHA256,
            "tables": {n: v for n, v in ext.K22_TABLE_SHA256.items() if v is not None}}


def reading_code() -> dict:
    """Which code made the reading: this script's and the shared module's sha256, the commit and its clean flag."""
    st = k18.stamp()
    return {"script": SCRIPT, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "outside_forms_sha256": hashlib.sha256((HERE / "outside_forms.py").read_bytes()).hexdigest(),
            "k21_sha256": hashlib.sha256((HERE / "k21_factor_forms.py").read_bytes()).hexdigest(),
            "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}


def reading(runs: Path, refs: Path = K18_TASKS, u3d: Path = k18.U3D_RECORDS) -> dict:
    """Apply the reading of row K22 to a run folder."""
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
    want_rec = {"dropped": dropped_arms(pre), "arms": list(arms), "ext": erec}
    rows: dict = {(w, n): [] for w, n in _cells()}
    for f in expected_files()[1:]:
        r, sha[f] = k21._read_json(runs / f)
        prov[f] = r["provenance"]
        world, rest = ("plain", f.split("/")[0][len("plain_"):]) if f.startswith("plain_") else ("dm", f[3:5])
        seed = int(f.rsplit("_", 1)[1][:-len(".json")])
        rec = r.get("pregate") or {}
        if rec.get("arms") != want_rec["arms"] or rec.get("dropped") != want_rec["dropped"]:
            why.append(f"{f}: arms or drops differ from pregate.json's")
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
                why.append(f"guard: {w} {n}: barred_form is off K18's committed barred_form" + (
                    f" by {d['mean_diff']:+.3f} in the mean and {d['max_abs_diff']:.3f} in one seed" if d else
                    " (a seed is missing)") + f" (limits {GUARD_MEAN} and {GUARD_SEED})")
    check = k18.dm_check({n: sorted(rows[("dm", n)], key=lambda q: q["seed"]) for n in TASKS}, u3d)
    tries = len(arms) - 1 + (0 if all(k18.guard_passes(d) for d in guard.values()) and len(guard) == 4 else 1)
    base = {"incomplete": False, "files_sha256": sha, "provenance": prov, "arms": list(arms), "pregate": pre,
            "ext": erec, "guard": guard, "dm_check_u3d_report_only": check, "tries_counted": tries}
    if why:
        return base | {"void": True, "why": why, "worlds": per, "outcome": "void: no file; T3's list is unchanged"}
    return (base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                    "worlds": per, "doubt": doubt_statement(per, arms)}
            | decide(per, arms, pre.get("corr_ln_all_ln_unrelated")))


def text(res: dict) -> str:
    lines = [f"Try K22: {res['outcome']}."]
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
                         f"{s['c']}, (d) {s['d']}); ahead in {s['ahead_cells'] or 'none'}; behind in "
                         f"{s['behind_cells'] or 'none'}; mean of the four paired means {s['mean_of_four']:+.2f}; "
                         f"gated (arm, base) {s['gated_vs_base']}")
        lines.append(f"  ranked: {res['ranked'] or 'none'}; file-2 rule {res['file2_rule']}; named: "
                     f"{res['files'] or 'none'}")
        d = res["doubt"]
        lines.append(f"  doubt statement (report only): {d['statement']}" + (
            "; rep_all - barred_form " + ", ".join(f"{k} {v:+.2f}" for k, v in d["means"].items()) if d["ran"]
            else ""))
    lines.append("guard (barred_form against K18's committed barred_form): " + "; ".join(
        f"{k} {k21._dist(d)}" for k, d in res["guard"].items()))
    lines.append("U3d's check arm in the dm world, reported only: " + "; ".join(
        f"{n} {k21._dist(c['distance'])}" for n, c in (res.get("dm_check_u3d_report_only") or {}).items()))
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    lines.append("g1 and g2 are two replicates of one knockout run; the controls and rep_all never qualify.")
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
                     f"{rc['outside_forms_sha256'][:12]}) at commit {rc['commit']}, clean {rc['tree_clean']}, "
                     f"machine {rc['machine']}")
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
