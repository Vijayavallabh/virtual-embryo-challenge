#!/usr/bin/env python
"""Try K21 (experiments/harness/TRIES.md; the row is placed before any run, and `run` and `pregate` refuse without it):
regularised forms of the barred factor read against the barred s1 on the two Gata4 replicates, in K18's two practice
worlds. A labelled reading that may name at most two T3 candidate files by the rule fixed in the row; not the adoption
rule.

Built on k18_built_files.py, imported unchanged: its stream (`t3_s2_form`), tasks (g1, g2), worlds (plain, dm), seeds,
the barred spec and factor (barred_spec, barred_factor, load_barred_for_plain, in_gene_order), the fold's factors
(fold_factors), the plain build (_build_plain), the depth-matched row (dm_row, _dm_task, quiet_selfcheck), the stamp,
the task-file writer, the folder rule, masking and the guard's distance. For a run, three of k18's module names are
replaced inside this process only (`k18_arms`): ARMS (this try's arms), build_arm (this try's arm builder, whose
barred_form branch is K18's own) and fold_factors (K18's, plus the permuted wild-type-only factors). So every arm draws
the barred s1's cells exactly as K18's barred_form did, and K18's barred_form is this try's base and guard.

The factor transforms are experiments/t3_decomp/factor_forms.py, the module the builder build_factor_form.py applies,
so a named file applies exactly the factor its practice arm read.

Arms (factor_forms.py has the formulas): barred_form (base, guard), top100, uni_top (control), cap90, uni_cap
(control), wtprod, wtprod_perm (control). The pre-gate (target-free, from the Mab21l2 table alone) drops a try with its
control when the form cannot be made or changes fewer than 10 genes; a dropped arm is not run and not counted.

  pregate --out DIR                                   the pre-gate: DIR/pregate.json and its lines (no score)
  run --world plain --task g1|g2 [--seeds 0-1] --out DIR   one plain task, seed by seed (resumable)
  run --world dm [--seeds ...] [--workers 8 --threads 2] --out DIR   both replicates, row by row (resumable)
  read --runs DIR --out DIR                           the fixed reading: reading.json and reading.txt

Layout of DIR: pregate.json, plain_g1/seed_<s>.json, plain_g2/seed_<s>.json, dm/setup.json, dm/<g1|g2>_seed_<s>.json,
each task file with its <file>.prov.json (commit, clean flag, machine, sha256), written by k18.write_task_file (by
rename, provenance first). A rerun into the same folder skips the seeds whose files are there; it is refused when the
folder holds a file of another commit or machine (k18.check_folder), and every task file's own stamp must be clean and
of the run's commit, else nothing is written for that seed. No log holds a score: the plain world hands the harness a
progress-only log; the depth-matched world prints K18's timing line per row, the frozen script's set-up lines and a
fixed self-check line; a crash prints its traceback with every decimal number masked (k18.masked).

The reading (fixed in the row; `reading`, `decide`):
  INCOMPLETE  pregate.json or any of the forty seed files missing: read opens nothing, writes nothing, exits 3.
  VOID        unless every file holds its seed and every arm that passed the pre-gate (the plain ones the floor too),
              all agree on the arms, come from a clean tree, one commit and one machine; and GUARD: barred_form within
              K17's tolerance (ten-seed mean 0.1, every seed 0.5) of K18's committed barred_form in all four cells
              (K18's task files at their pinned sha256). Expected distance 0 on the same machine; nonzero seeds listed.
  AHEAD / BEHIND in a cell: arm minus barred_form, paired by seed, mean > 0.3 with >= 8 of 10 positive / mean < -0.3.
  QUALIFIED   (top100, cap90, wtprod only): (a) AHEAD in a plain cell and in a dm cell, and on g1 in some world and on
              g2 in some world; (b) BEHIND in none; (c) in every cell its gated seeds at most barred_form's plus 1;
              (d) in every cell where it is AHEAD, it beats its control by a mean > 0.3 with >= 8 of 10 positive.
  NAMES       rank the qualified by the number of AHEAD cells, then the mean of the four paired means (ties: top100,
              cap90, wtprod); file 1 the first form on base B, file 2 the second (decision log 2026-11-20 (q)'s B).
              None if none qualifies or the run is void.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
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
import factor_forms as ff  # noqa: E402  (the transforms, shared with build_factor_form.py)
import k7_wt_revert as k7  # noqa: E402  (the permutation seed, the exponent, _paired, MIN_GAIN and MIN_AHEAD)
import k8_s2_form as k8  # noqa: E402  (the upload's central composition, for the record)
import k16_revert_dose as k16  # noqa: E402  (move_by_group: K8b's per-group step)
import k18_built_files as k18  # noqa: E402  (the worlds, tasks, stream, seeds, build path and scoring; unchanged)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.evaluate import SEVERITY_GATED  # noqa: E402

ROW_SCRIPT = "`experiments/t3_decomp/k21_factor_forms.py`"  # the ledger row names this script
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"
STREAM = k18.STREAM  # "t3_s2_form"
BASE = k18.BARRED  # "barred_form"
ARMS = (BASE, "top100", "uni_top", "cap90", "uni_cap", "wtprod", "wtprod_perm")
VECTOR_OF = {BASE: "identity", "top100": "top100", "uni_top": "uni_top", "cap90": "cap90", "uni_cap": "uni_cap"}
NAMEABLE = ff.TRIES  # ("top100", "cap90", "wtprod"), also the tie order
CONTROL_OF = dict(ff.CONTROL_OF)
TASKS = k18.TASKS
WORLDS = k18.WORLDS
MIN_GAIN, MIN_AHEAD = k7.MIN_GAIN, k7.MIN_AHEAD  # 0.3 and 8 of 10
MAX_EXTRA_GATED = 1
GUARD_MEAN, GUARD_SEED = k18.GUARD_MEAN, k18.GUARD_SEED  # 0.1 and 0.5
K18_TASKS = HERE / "results" / "k18_built_files" / "tasks"
K18_SHA256 = {"plain_g1.json": "5ec09ef583862ab47ab078e72c7aed807c89b4d9c234c57e3b0727ff9e460ed9",
              "plain_g2.json": "cfc5ede7e7185ba5e16d4774d975b9e140e60bbd29736fb85c0fdb8e99f4634f",
              "dm.json": "1ef3e26347c65f915c1630130da02a655c155b912392180001cba315fe135b60"}
RESULTS = HERE / "results" / "k21_factor_forms"


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row


def row_present(path: Path | None = None) -> str:
    """The id of the TRIES.md row that names this script (a table line starting '| <id> ' and holding ROW_SCRIPT);
    refuses (SystemExit) when there is none or more than one."""
    path = TRIES_MD if path is None else path
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.startswith("| ") and ROW_SCRIPT in ln]
    ids = sorted({ln[2:].split(" ", 1)[0] for ln in rows})
    if len(ids) != 1:
        sys.exit(f"refused: {'no' if not ids else 'more than one'} row of {Path(path).name} names {ROW_SCRIPT} "
                 f"({ids}); the row is committed before any run")
    return ids[0]


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K18's fold_factors (the original, called unchanged) plus the wild-type-only factors permuted as
    k16.permuted_factors does it: within deciles of `basis`, seed k7.PERM_SEED, one permutation for every group and
    the pooled factor."""
    F = _K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    F["per_perm"] = k18.permuted(per, basis)
    F["phi_perm"] = harness.t3_decile_permutation(phi, basis, k7.PERM_SEED)
    return F


def arm_factor(arm: str, phi_m: np.ndarray) -> np.ndarray:
    """A vector arm's factor from phi_M in the reference's gene order (barred_form: phi_M itself, the same array)."""
    phi, _ = ff.vector_form(VECTOR_OF[arm], phi_m)
    return phi


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """One arm's drawn cells. Every arm draws from stable_rng(seed, STREAM) through bs.construct with the barred
    spec's keywords (prediction_inputs: the central composition, w = 1, no depth ratio or stage mixing), so every arm
    holds the barred s1's rows. A vector arm passes its factor in place of phi (barred_form: K18's arm exactly). wtprod
    and its control draw with no factor, then move each cell by group through k16.move_by_group at k7.EXPONENT on
    stable_rng(seed, STREAM:grouped) with per'_g = phi_M**2 * per_g (fallback phi_M**2 * pooled)."""
    if arm not in ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    if arm in VECTOR_OF:
        kw["phi"] = arm_factor(arm, F["phi_barred"])
        return bs.construct(ref, n, rng, **kw)
    kw["phi"] = None
    built = bs.construct(ref, n, rng, **kw)
    per, pooled = (F["per"], F["phi"]) if arm == "wtprod" else (F["per_perm"], F["phi_perm"])
    per2, fallback2 = ff.wtprod_factors(F["phi_barred"], per, pooled)
    return k16.move_by_group(bs, ref, built, F["group"], per2, fallback2, k7.EXPONENT,
                             bs.t3d.stable_rng(seed, f"{STREAM}:grouped"))


_K18_ORIGINAL = {"ARMS": k18.ARMS, "build_arm": k18.build_arm, "fold_factors": k18.fold_factors}


@contextmanager
def k18_arms(arms: tuple[str, ...]):
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's, restored on exit. K18's
    file is not changed; its dm_row and _build_plain then build this try's arms with K18's own set-up and scoring."""
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
                                         "factor": "the barred s1's (build_submission.py --branch barred --submission "
                                                   f"s1) in the form {arm} (experiments/t3_decomp/factor_forms.py)",
                                         "perm_seed": k7.PERM_SEED if arm == "wtprod_perm" else None},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


# ---------------------------------------------------------------------------------------------------------------------
# Running


def parse_seeds(text_: str | None) -> list[int]:
    """'0-1' or '3' or '0,4,7' or None (all ten); only seeds of vec.adoption.SEEDS."""
    if text_ is None:
        return list(SEEDS)
    out: list[int] = []
    for part in str(text_).split(","):
        a, _, b = part.partition("-")
        lo, hi = int(a), int(b) if b else int(a)
        out += list(range(lo, hi + 1))
    if not out or any(s not in SEEDS for s in out) or len(set(out)) != len(out):
        raise SystemExit(f"--seeds {text_!r}: seeds must be distinct members of {list(SEEDS)}")
    return sorted(out)


def progress_log(world: str, task: str):
    """The harness's per-seed line lists every arm's score; this prints the world, the task and the seed only."""

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"K21 {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


def _stamp_ok(st: dict) -> None:
    """A task file is written only if the tree is still clean and at the run's commit."""
    now = k18.stamp()
    if now["tree_clean"] is not True or now["commit"] != st["commit"]:
        raise _Unclean(f"the tree is now commit {now['commit']}, clean {now['tree_clean']}, against {st['commit']} "
                       "clean at the start: no task file written (INCOMPLETE; rerun on a clean tree at the commit)")


class _Unclean(RuntimeError):
    pass


def _pregate_record(pre: dict, arms: tuple[str, ...]) -> dict:
    return {"dropped": pre["dropped"], "arms": list(arms)}


def run_plain(a, st: dict, out: Path) -> None:
    """One plain task, seed by seed: each seed through harness.run_task (K8's harness) with this try's candidates,
    written as plain_<task>/seed_<s>.json; seeds already written at this commit on this machine are skipped."""
    folder = out / f"plain_{a.task}"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    pre = ff.pregate(k18._BARRED["phi"])
    arms = ff.active_arms(pre, ARMS)
    with k18_arms(arms):
        for s in parse_seeds(a.seeds):
            path = folder / f"seed_{s}.json"
            if path.exists():
                print(f"K21 plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            res = harness.run_task(TASKS[a.task], candidates(arms), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log("plain", a.task))
            got = res["provenance"]
            if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
                raise _Unclean(f"the harness's stamp says commit {got.get('commit')}, clean {got.get('tree_clean')}: "
                               "no task file written (INCOMPLETE)")
            _stamp_ok(st)
            res["provenance"]["machine"] = st["machine"]
            res |= {"world": "plain", "task_source": "experiments/t3_decomp/k21_factor_forms.py", "seed": s,
                    "pregate": _pregate_record(pre, arms), "factor_record": k7._PREP["fold"]["record"],
                    "barred_factor": rec}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)


def run_dm(a, st: dict, out: Path) -> None:
    """The depth-matched world, set up exactly as k18.run_dm sets it up (the same functions in the same order), then
    each (replicate, seed) row through k18.dm_row with this try's arms, written as dm/<task>_seed_<s>.json as it
    finishes; rows already written at this commit on this machine are skipped."""
    import multiprocessing as mp

    import orth_diag as od

    folder = out / "dm"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    seeds = parse_seeds(a.seeds)
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
        pre = ff.pregate(phi_b)
        arms = ff.active_arms(pre, ARMS)
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
                sys.exit(f"{sp} records another set-up (inputs, factor or pre-gate): not run")
        else:
            k18.write_task_file(sp, setup | {"provenance": st})
        with k18_arms(arms):
            F = k18.fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b,
                                 spec)
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
                                                 "task_source": "experiments/t3_decomp/k21_factor_forms.py"})
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


DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing written. Row K21: commit or restore first; INCOMPLETE, not VOID")


def cmd_run(a) -> None:
    if a.world == "plain" and a.task not in TASKS:
        sys.exit("--world plain needs --task g1 or g2")
    row_present()
    st = k18.stamp()
    if not st["tree_clean"]:  # before anything loads
        sys.exit(DIRTY)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        (run_plain if a.world == "plain" else run_dm)(a, st, out)
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seeds written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)


def pregate_lines(pre: dict) -> list[str]:
    def fmt(v):
        return f"{v:.6g}" if isinstance(v, float) else str(v)

    lines = [f"pre-gate (the Mab21l2 table alone; no score): {pre['n_genes']} genes, "
             f"{pre['genes_at_clip_phi_m']} at the clip; phi_M " + ", ".join(
                 f"{k} {fmt(v)}" for k, v in pre["phi_m"].items())]
    for arm, r in pre["arms"].items():
        lines.append(f"  {arm}: " + ", ".join(f"{k} {fmt(v)}" for k, v in r.items() if not isinstance(v, dict)))
    lines.append("  dropped: " + ("; ".join(f"{k} ({v})" for k, v in pre["dropped"].items()) or "none"))
    return lines


def cmd_pregate(a) -> None:
    """The pre-gate before any scoring: the barred factor through the builder's loader (WT E8.75's gene order) and the
    wild-type-only factors through the builder's wt_revert_factors (for wtprod's printed range); DIR/pregate.json."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        sys.exit(DIRTY)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "pregate.json"
    if path.exists():
        sys.exit(f"{path} exists")
    try:
        bs = harness._t3_head()
        rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
        phi_m, genes = k18._BARRED["phi"], k18._BARRED["genes"]
        wt = Path(a.data) / "T2" / "E8.75.h5ad"
        series = bs.wt_revert_series(wt, {n: Path(a.data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                          if n != wt.stem})
        pooled, per, _, wrec = bs.wt_revert_factors(wt, series, genes)
        pre = ff.pregate(phi_m, per, pooled)
        res = {"pregate": pre, "arms": list(ff.active_arms(pre, ARMS)), "row": rid, "barred_factor": rec,
               "wt_revert": {k: v for k, v in wrec.items() if k != "series"},
               "wt_series_sha256": {n: v["sha256"] for n, v in wrec["series"].items()}, "provenance": st}
        k18.write_task_file(path, res)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print("\n".join(pregate_lines(pre)), flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _gated(sev: list[dict]) -> int:
    return int(sum(v["state"] in SEVERITY_GATED for v in sev))


def summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """Per arm: the score's mean and sd, each metric's mean skill, the severity states, seeds scored and gated, the
    median beta and r2. Paired (k7._paired): each arm against barred_form, each try against its control, barred_form
    against the floor (plain). Per-term: each arm's mean per-seed skill difference against barred_form."""
    out: dict = {"arms": {}, "paired": {}, "per_term": {}}
    for x in arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "severity_states": {v: sum(z["state"] == v for z in sev) for v in sorted({z["state"] for z in sev})},
            "severity_scored": int(sum(v["state"] == "scored" for v in sev)), "severity_gated": _gated(sev),
            "median_beta": k18._median([v.get("beta") for v in sev]),
            "median_r2": k18._median([v.get("r2") for v in sev])}
    pairs = [(x, BASE) for x in arms if x not in (BASE, harness.FLOOR)]
    pairs += [(t, c) for t, c in CONTROL_OF.items() if t in arms and c in arms]
    if harness.FLOOR in arms:
        pairs.append((BASE, harness.FLOOR))
    out["paired"] = {f"{x} - {y}": k7._paired(rows, x, y) for x, y in pairs}
    for x in arms:
        if x in (BASE, harness.FLOOR):
            continue
        ms = rows[0]["skills"][BASE]
        out["per_term"][x] = {m: float(np.mean([q["skills"][x][m] - q["skills"][BASE][m] for q in rows])) for m in ms}
    return out


def _cells():
    return [(w, n) for w in WORLDS for n in TASKS]


def decide(per: dict, arms: tuple[str, ...]) -> dict:
    """The rule fixed in the row, on summarise()'s output per (world, task)."""
    status: dict = {}
    for t in NAMEABLE:
        if t not in arms:
            status[t] = {"ran": False, "qualified": False}
            continue
        c = CONTROL_OF[t]
        d = {(w, n): per[w][n]["paired"][f"{t} - {BASE}"] for w, n in _cells()}
        ahead = {k: bool(v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD) for k, v in d.items()}
        behind = {k: bool(v["mean"] < -MIN_GAIN) for k, v in d.items()}
        gated = {k: (per[k[0]][k[1]]["arms"][t]["severity_gated"], per[k[0]][k[1]]["arms"][BASE]["severity_gated"])
                 for k in d}
        ctrl = {k: per[k[0]][k[1]]["paired"][f"{t} - {c}"] for k in d}
        a = (any(ahead[("plain", n)] for n in TASKS) and any(ahead[("dm", n)] for n in TASKS)
             and any(ahead[(w, "g1")] for w in WORLDS) and any(ahead[(w, "g2")] for w in WORLDS))
        b = not any(behind.values())
        cc = all(g <= g0 + MAX_EXTRA_GATED for g, g0 in gated.values())
        dd = all(ctrl[k]["mean"] > MIN_GAIN and ctrl[k]["ahead"] >= MIN_AHEAD for k in d if ahead[k])
        status[t] = {"ran": True, "qualified": bool(a and b and cc and dd), "a": bool(a), "b": bool(b), "c": bool(cc),
                     "d": bool(dd), "ahead_cells": [f"{w} {n}" for (w, n), v in ahead.items() if v],
                     "behind_cells": [f"{w} {n}" for (w, n), v in behind.items() if v],
                     "mean_of_four": float(np.mean([v["mean"] for v in d.values()])),
                     "gated_vs_base": {f"{w} {n}": list(v) for (w, n), v in gated.items()}}
    q = [t for t in NAMEABLE if status[t]["qualified"]]
    ranked = sorted(q, key=lambda t: (-len(status[t]["ahead_cells"]), -status[t]["mean_of_four"], NAMEABLE.index(t)))
    names = ranked[:2]
    return {"status": status, "ranked": ranked, "names": names,
            "files": [f"{t}-B (base B of decision log 2026-11-20 (q); build_factor_form.py --form {t}-<B>)"
                      for t in names],
            "outcome": ("names " + " and ".join(f"file {i + 1}: {t} on base B" for i, t in enumerate(names))
                        + "; each after the C list, under the R14/R17 conditions" if names else
                        "no file: nothing qualifies; T3's list is unchanged")}


def _read_json(p: Path) -> tuple[dict, str]:
    raw = p.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def expected_files() -> list[str]:
    return (["pregate.json"] + [f"plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def k18_rows(refs: Path = K18_TASKS) -> tuple[dict, list[str]]:
    """K18's committed barred_form rows per cell, from its task files at their pinned sha256."""
    why, out = [], {}
    for f, sha in K18_SHA256.items():
        p = refs / f
        if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest() != sha:
            why.append(f"K18's committed {f} is missing or not at its pinned sha256")
            continue
        r = json.loads(p.read_text())
        if f == "dm.json":
            for n in TASKS:
                out[("dm", n)] = r["rows"].get(n, [])
        else:
            out[("plain", f[len("plain_"):-len(".json")])] = r["rows"]
    return out, why


INCOMPLETE = ("not a reading: pregate.json or a seed's task file is missing (a crash, a refusal such as a dirty tree, "
              "out of memory or a kill); nothing was read. Row K21: if the cause is outside the code, rerun the "
              "missing seeds into this folder at the run's commit on the same machine (the runs skip what is there); "
              "after a committed code fix, rerun everything at the new commit into a fresh folder and never read this "
              "one")


def reading(runs: Path, refs: Path = K18_TASKS, u3d: Path = k18.U3D_RECORDS) -> dict:
    """Apply the reading of row K21 to a run folder."""
    runs = Path(runs)
    missing = [f for f in expected_files() if not (runs / f).exists()]
    if missing:
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    why: list[str] = []
    sha, prov = {}, {}
    pre_file, sha["pregate.json"] = _read_json(runs / "pregate.json")
    prov["pregate.json"] = pre_file["provenance"]
    pre = pre_file["pregate"]
    arms = ff.active_arms(pre, ARMS)
    if list(arms) != pre_file["arms"]:
        why.append("pregate.json: its arms are not the pre-gate's")
    rows: dict = {(w, n): [] for w, n in _cells()}
    for f in expected_files()[1:]:
        r, sha[f] = _read_json(runs / f)
        prov[f] = r["provenance"]
        world, rest = ("plain", f.split("/")[0][len("plain_"):]) if f.startswith("plain_") else ("dm", f[3:5])
        seed = int(f.rsplit("_", 1)[1][:-len(".json")])
        rec = r.get("pregate") or {}
        if rec.get("arms") != list(arms) or rec.get("dropped") != pre["dropped"]:
            why.append(f"{f}: arms or drops differ from pregate.json's")
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
    ref_rows, kwhy = k18_rows(refs)
    why += kwhy
    per: dict = {"plain": {}, "dm": {}}
    guard: dict = {}
    for (w, n), rr in rows.items():
        rr = sorted(rr, key=lambda q: q["seed"])
        if [q["seed"] for q in rr] != list(SEEDS):
            why.append(f"{w} {n}: not the ten seeds")
            continue
        per[w][n] = summarise(rr, (harness.FLOOR, *arms) if w == "plain" else arms)
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
            "guard": guard, "dm_check_u3d_report_only": check, "tries_counted": tries}
    if why:
        return base | {"void": True, "why": why, "worlds": per,
                       "outcome": "void: no file; T3's list is unchanged"}
    return base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                   "worlds": per} | decide(per, arms)


def _f(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def _dist(d: dict | None) -> str:
    if d is None:
        return "n/a"
    return (f"mean {d['mean_diff']:+.4f}, largest seed {d['max_abs_diff']:.4f}, seeds differing "
            f"{d['seeds_off'] or 'none'}")


def text(res: dict) -> str:
    lines = [f"Try K21: {res['outcome']}."]
    if res.get("incomplete"):
        miss = res["missing"]
        return "\n".join([*lines, f"missing: {len(miss)} file(s), e.g. " + ", ".join(miss[:6])]) + "\n"
    lines.append(f"arms run: {', '.join(res['arms'])}; tries counted: {res['tries_counted']}")
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
        lines.append(f"  ranked: {res['ranked'] or 'none'}; named: {res['files'] or 'none'}")
    lines.append("guard (barred_form against K18's committed barred_form): " + "; ".join(
        f"{k} {_dist(d)}" for k, d in res["guard"].items()))
    lines.append("U3d's check arm in the dm world, reported only: " + "; ".join(
        f"{n} {_dist(c['distance'])}" for n, c in (res.get("dm_check_u3d_report_only") or {}).items()))
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    lines.append("g1 and g2 are two replicates of one knockout run; the controls never qualify.")
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
                             f"{_f(s['median_beta'])}, r2 {_f(s['median_r2'])}; skills "
                             + ", ".join(f"{m} {x:.3f}" for m, x in s["skills"].items()))
            for k, d in v["paired"].items():
                lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
            for a, terms in v["per_term"].items():
                lines.append(f"    per-term skill, {a} - {BASE}: " + ", ".join(f"{m} {x:+.3f}" for m, x in
                                                                              terms.items()))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    try:
        res = reading(Path(a.runs))
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
    except Exception:  # noqa: BLE001  (masked traceback, nothing written)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if res.get("incomplete"):
        print(t, end="")
        sys.exit(3)
    out = Path(a.out)
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
