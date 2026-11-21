#!/usr/bin/env python
"""Row E53 of experiments/harness/TRIES.md (as specified by a score-blind judge as moonshot M2, "Geodesic overlay
(LDDMM midpoint of the two brackets)"; the specification names its row E51, renamed E53 because labels collided; placed
by the main session before any run): the kill test of the geodesic overlay on the embryo's exact-position file, on the
nine practice tasks. Two tries, both J5 diagnostics: E53a the arm, E53b its rotated-start placebo.

Arms (harness candidates, vec.harness.GEODESIC_CANDIDATES; each is fgw_shape0x_ov's rows, order, count, expression and
pool weights with positions on one LDDMM geodesic, vec.shape_interp.geodesic_overlay and vec.lddmm):
  E53a fgw_shape0x_geo      the arm: the later bracket in the earlier frame by the build's own FGW rotation R
  E53b fgw_shape0x_geo_rot  the placebo: a uniform proper rotation Q from stream [build seed, 50] in place of R
Rescored, not counted: C0 = fgw_shape0x and E49c = fgw_shape0x_ov.

Path. harness.run_task on each task of nine_tasks.all_tasks(), seeds 0-9, bake protocol, build seed BUILD_SEED_OFFSET +
seed, under build "main" (primary) and "prefix" (reported, deciding nothing). The arm, the placebo and E49c are wrapped
as E49 wraps its forms (e49_shape_interp_read.guarded): a refusal (vec.shape_interp.ShapeRefused or the harness's cloud
guard) is recorded and C0's file of that seed stands in; any other error is raised.

Steps (each refuses unless row E53 is in TRIES.md, the tree is clean and HEAD is the run commit; none meets a test
board):
  describe --pair NAME --commit SHA  target-free, not a try: the arm's fit at build seeds 1,000,001 to 1,000,010 on one
                                     pair of brackets ("test": E7.5 and E8.0 with the validation board's panel and count
                                     rule at t = 1/2; or a task's two brackets with its panel board, at its t). Reads no
                                     target. One file per pair.
  pregate --commit SHA               the gate clauses g1 to g3 on the test pair's ten draws (the nine tasks' values are
                                     reported); exit 0 pass, 4 fail (nothing runs, nothing counts)
  j7 --commit SHA                    J7's dry run: the 0.01 reproduction tolerance against E48's committed rescore
                                     (experiments/coords/results/e48_row_order/gate, its worst deviations from E47);
                                     the result is quoted in the reading, deciding nothing
  gate --task NAME --commit SHA      the reproduction gate: the floor and C0 under main against E47's committed per-seed
                                     scores (0.01 per seed and per ten-seed mean, as E49's gate), and E49c against E49's
                                     committed reading (its per-task ten-seed mean and SD of E49c minus C0, 0.01 each);
                                     exit 3 on a miss
  run --task NAME --commit SHA       C0, E49c, the arm and the placebo, both builds, seed by seed (resumable); refused
                                     unless the pregate, j7 and all nine gate files pass at the same commit
  read --commit SHA                  the reading into experiments/coords/results/e53_geodesic/ (reading.json,
                                     reading.txt), never overwriting
The script prints only step lines ("describe <pair> done", "pregate pass|fail", "j7 within|beyond", "gate <task>
pass|fail", "<task> seed <s> done", "wrote <path>"), never a count. It reads no test board, test file or test anchor.

Target-free gate, fixed in the row (arm only, before any score): on the test pair,
  g1  the data term at t = 1 falls to at most 0.5 of its value at the rigid start in at least 8 of 10 draws;
  g2  the smallest Jacobian determinant of the map over the drawn rows at t is above 0.05 in every draw;
  g3  the moved union's thinnest normalised SD lies between the two brackets' own thinnest SDs in at least 8 of 10
      draws.
If a clause fails, nothing runs and nothing counts.

Pass line, fixed in the row (build main; merged over five truths at equal weight as E49: E8.0 = half, late; E8.25_late
= new1, x1; E8.5 = new2, x2, x3; E7.5 = e75; E7.25 = e725; test-like = half, late, new1, new2). All must hold:
  (a) merged board gain of the arm over E49c at least +1.0, the truth means positive on at least 4 of 5 truths, and the
      merged gain at least +0.5 with any single truth left out;
  (b) test-like board gain over E49c at least +1.0;
  (c) gated shape_growth group skill gain at least +0.04 over C0 and at least +0.02 over E49c;
  (d) local_organisation skill not below E49c's by more than 0.01 merged;
  (e) arm minus placebo on the ungated shape S = (d2_shape skill + occupancy_dice skill) / 2 at least +0.03 merged,
      with at least 8 of 10 seeds ahead;
  (f) the arm refused (cloud guard or form) on at most 5 of 90 task-seeds.
Outcome: if (a) to (f) all hold, the candidate file is fixed in the decision log before any score of it and goes to the
captain as a priced question beside E49c; otherwise no file.
Void (read as nothing): a missing task, seed, build or arm; a result from another commit or an unclean tree; a
pregate, j7 or gate file missing or failing; any arm whose count or expression matrix (in row order) differs from C0's;
an expression-only skill off C0's by more than 1e-9; scale_log_ratio's skill off C0's by more than 1e-6 where neither
file was jittered; a file or a refusal that differs between the builds; a test board (E49's void list).
Reported, deciding nothing: the per-truth table; e75 and e725 separately; d2_shape and occupancy_dice per task; the gate
values; the prefix build; the placebo's and E49c's refusals.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))
sys.path.insert(0, str(REPO / "experiments" / "harness"))
sys.path.insert(0, str(REPO / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
import e47_floor_exact as e47  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import evaluate as ev  # noqa: E402
from vec import harness, lddmm  # noqa: E402
from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402


def _load_e49():
    """E49's script, for its wrappers, series, merge and E47 gate (imported, not copied)."""
    name = "e49_shape_interp_read"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, HERE / "e49_shape_interp_read.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


e49 = _load_e49()

ROW_IDS = ("E53a", "E53")  # the row's first word: "E53a and E53b" as drafted, or "E53"
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
C0, E49C = "fgw_shape0x", "fgw_shape0x_ov"
ARM, PLACEBO = "fgw_shape0x_geo", "fgw_shape0x_geo_rot"
TRY_IDS = {ARM: "E53a", PLACEBO: "E53b"}
FORM_OF = {**harness.SHAPE_CANDIDATES, **harness.GEODESIC_CANDIDATES}
WRAPPED = (E49C, ARM, PLACEBO)  # guarded: a refusal is recorded and C0's file stands in
NAMES = (C0, *WRAPPED)
BUILDS = ("main", "prefix")
TEST_LIKE = ("half", "late", "new1", "new2")
SINGLE_TRUTH_TASKS = ("e75", "e725")
GROUP_SHAPE, GROUP_LOCAL = "shape_growth", "local_organisation"
SHAPE_TERMS = e49.SHAPE_TERMS
EXPR_METRICS = e49.EXPR_METRICS
EXPR_TOL, SCALE_TOL = e49.EXPR_TOL, e49.SCALE_TOL
# the pass line
A_MERGED, A_TRUTHS, A_LEAVE_ONE_OUT = 1.0, 4, 0.5
B_TEST_LIKE = 1.0
C_OVER_C0, C_OVER_E49C = 0.04, 0.02
D_LOCAL = 0.01
E_S, E_AHEAD = 0.03, 8
F_MAX_REFUSED = 5
# the target-free gate
G1_RATIO, G1_DRAWS = 0.5, 8
G2_JAC = 0.05
G3_DRAWS = 8
DESCRIBE_SEEDS = tuple(harness.BUILD_SEED_OFFSET + 1 + i for i in range(10))
DESCRIBE_BOARD = "T2:embryo:val_interp"
TEST_PAIR = {"before": "T2/E7.5.h5ad", "after": "T2/E8.0.h5ad"}
TEST_T = 0.5
PAIRS = ("test", *nt.all_tasks())
# the reproduction gate
REPRO_TOL = 0.01  # board points (J7)
E47_DIR, E47_RUN_COMMIT, E47_FILES_COMMIT = e49.E47_DIR, e49.E47_RUN_COMMIT, e49.E47_FILES_COMMIT
E49_READING = REPO / "experiments" / "coords" / "results" / "e49_shape_interp" / "reading.json"
E49_FILES_COMMIT = "53ba6cc"
E48_GATE_DIR = REPO / "experiments" / "coords" / "results" / "e48_row_order" / "gate"
E48_FILES_COMMIT = "af83195"
BASES = ("c902b90",)  # main when this branch was cut
SCRIPT = "experiments/coords/e53_geodesic_read.py"
RUNS = REPO / "runs" / "e53"
RESULTS = REPO / "experiments" / "coords" / "results" / "e53_geodesic"


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def has_row(path: Path | None = None) -> bool:
    rows = n12.tries_rows(TRIES if path is None else path)
    return any(r in rows for r in ROW_IDS)


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E53 at `commit` (empty: it may)."""
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    for base in BASES:
        if _git("merge-base", "--is-ancestor", base, "HEAD").returncode != 0:
            why.append(f"HEAD does not contain {base}")
    for d, c in ((E47_DIR, E47_FILES_COMMIT), (E49_READING, E49_FILES_COMMIT), (E48_GATE_DIR, E48_FILES_COMMIT)):
        if _git("diff", "--quiet", c, "HEAD", "--", str(d.relative_to(REPO))).returncode != 0:
            why.append(f"{d.relative_to(REPO)} differs from {c}'s")
    return why


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (ev, harness, tb, si, lddmm):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    if not has_row():
        sys.exit("E53 has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E53 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E53 does not run on one")
    return t


def _write(path: Path, rec: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec, indent=1, default=float))
    os.replace(tmp, path)


def _read_at(path: Path, commit: str) -> tuple[dict | None, list[str]]:
    if not path.exists():
        return None, [f"no {path.name}"]
    d = json.loads(path.read_text())
    return d, ([] if d.get("commit") == commit else [f"{path.name} is from another commit"])


# ------------------------------------------------------------------------------------------------ target-free gate
def task_fraction(t) -> float:
    """A task's calendar fraction, as harness._fgw computes it (from the stage days of its files)."""
    lo, hi = harness.STAGE_DAYS[t.inputs["before"]], harness.STAGE_DAYS[t.inputs["after"]]
    return (harness.STAGE_DAYS[t.target] - lo) / (hi - lo)


def gate_values(record: dict) -> dict:
    """The arm's gate quantities from one interpolate_clouds record (flat values; no count)."""
    return {"data_ratio": float(record["shape_interp_lddmm_data_ratio"]),
            "jac_min": float(record["shape_interp_jac_min"]),
            "thin_sd": float(record["shape_interp_thin_sd"]),
            "thin_sd_before": float(record["shape_interp_thin_sd_before"]),
            "thin_sd_after": float(record["shape_interp_thin_sd_after"]),
            "thin_between": bool(record["shape_interp_thin_between"]),
            "lbfgs_n_iter": int(record.get("shape_interp_lddmm_lbfgs_n_iter", -1)),
            "p0_norm_mean": float(record.get("shape_interp_lddmm_p0_norm_mean", float("nan")))}


def gate_clauses(draws: list[dict]) -> dict:
    """g1 to g3 on a pair's draws (each a gate_values dict, or {"refused": why})."""
    ok = [d for d in draws if "refused" not in d]
    g1 = sum(d["data_ratio"] <= G1_RATIO for d in ok)
    g2 = all(d["jac_min"] > G2_JAC for d in ok) and len(ok) == len(draws)
    g3 = sum(bool(min(d["thin_sd_before"], d["thin_sd_after"]) <= d["thin_sd"]
                  <= max(d["thin_sd_before"], d["thin_sd_after"])) for d in ok)
    c = {"g1": bool(g1 >= G1_DRAWS), "g2": bool(g2), "g3": bool(g3 >= G3_DRAWS)}
    return {"clauses": c, "passes": all(c.values()), "g1_draws_at_or_below": int(g1), "g3_draws_between": int(g3),
            "refused_draws": len(draws) - len(ok), "jac_min": float(min((d["jac_min"] for d in ok), default=np.nan)),
            "data_ratio_max": float(max((d["data_ratio"] for d in ok), default=np.nan))}


def _pair_inputs(name: str, data: Path):
    from vec.boards import load_boards
    from vec.io import COORD_KEY, dense
    from vec.prep import load_stage

    boards = load_boards()
    if name == "test":
        board, files, t = boards[DESCRIBE_BOARD], TEST_PAIR, TEST_T
    else:
        task = _task(name)
        files = {r: task.inputs[r] for r in ("before", "after")}
        board, t = boards[task.panel_board], task_fraction(task)
    if e49._is_test_board(board):
        sys.exit("describe reads a validation board's panel and count rule only")
    genes = list(board.genes)
    st = {r: load_stage(data / f, genes=genes, drop_duplicates=True) for r, f in files.items()}
    C = [np.asarray(st[r].obsm[COORD_KEY], float)[:, :3] for r in ("before", "after")]
    E = [dense(st[r][:, genes].X).astype(np.float32) for r in ("before", "after")]
    return C, E, t, harness.count_for(board, st["before"].n_obs), board.key, files


def describe_pair(C, E, t: float, n: int, seeds=DESCRIBE_SEEDS, form: str = "geodesic_overlay") -> list[dict]:
    """The arm's gate values at each seed on one pair of brackets (no target)."""
    out = []
    for seed in seeds:
        try:
            res = tb.interpolate_clouds(C[0], C[1], E[0], E[1], t, n, seed, target_rms=1.0, t_pos=0.0,
                                        positions=tb.POSITIONS[1], shape_form=form)
            out.append({"seed": int(seed), **gate_values(res.record)})
        except si.ShapeRefused as e:
            out.append({"seed": int(seed), "refused": str(e)})
    return out


def cmd_describe(a) -> None:
    _refuse_unless_ready(a.commit)
    if a.pair not in PAIRS:
        sys.exit(f"unknown pair {a.pair}; known: {PAIRS}")
    out = Path(a.runs) / "describe" / f"{a.pair}.json"
    if out.exists():
        sys.exit(f"{out} exists")
    C, E, t, n, board, files = _pair_inputs(a.pair, Path(a.data))
    rec = {"commit": a.commit, "script": SCRIPT, "pair": a.pair, "files": files, "panel_and_count_rule": board,
           "t": float(t), "seeds": list(DESCRIBE_SEEDS), "draws": describe_pair(C, E, t, n)}
    _write(out, rec)
    print(f"describe {a.pair} done", flush=True)


def pregate(runs: Path, commit: str) -> dict:
    """g1 to g3 on the test pair decide; every pair's values are kept (the nine tasks reported only)."""
    why, per = [], {}
    for p in PAIRS:
        d, w = _read_at(Path(runs) / "describe" / f"{p}.json", commit)
        why += [f"describe {p}: {x}" for x in w]
        if d is not None and not w:
            if [x["seed"] for x in d["draws"]] != list(DESCRIBE_SEEDS):
                why.append(f"describe {p}: seeds are not the gate's")
                continue
            per[p] = {"t": d["t"], "draws": d["draws"], **gate_clauses(d["draws"])}
    passes = not why and per.get("test", {}).get("passes", False)
    return {"commit": commit, "why": why, "pairs": per, "decides": "test", "pass": bool(passes)}


def cmd_pregate(a) -> None:
    _refuse_unless_ready(a.commit)
    out = Path(a.runs) / "pregate.json"
    if out.exists():
        sys.exit(f"{out} exists")
    rec = pregate(Path(a.runs), a.commit)
    _write(out, rec)
    print(f"pregate {'pass' if rec['pass'] else 'fail'}", flush=True)
    sys.exit(0 if rec["pass"] else 4)


# ------------------------------------------------------------------------------------------------ J7 dry run
def j7(gate_dir: Path = E48_GATE_DIR, tol: float = REPRO_TOL) -> dict:
    """The worst deviations E48's committed rescore found against E47 (per seed, per ten-seed mean), against tol."""
    worst = {"mean": 0.0, "seed": 0.0}
    files = sorted(Path(gate_dir).glob("*.json"))
    why = [] if len(files) == len(nt.all_tasks()) else ["E48's gate files are not one per task"]
    for p in files:
        w = json.loads(p.read_text()).get("worst") or {}
        if w.get("mean") is None or w.get("seed") is None:
            why.append(f"{p.name}: no worst deviation")
            continue
        worst = {k: max(worst[k], float(w[k])) for k in worst}
    within = not why and worst["mean"] <= tol and worst["seed"] <= tol
    return {"source": str(Path(gate_dir).relative_to(REPO)) if Path(gate_dir).is_relative_to(REPO) else str(gate_dir),
            "files_commit": E48_FILES_COMMIT, "tolerance": tol, "worst": worst, "why": why, "within": bool(within)}


def cmd_j7(a) -> None:
    _refuse_unless_ready(a.commit)
    out = Path(a.runs) / "j7.json"
    if out.exists():
        sys.exit(f"{out} exists")
    rec = {**j7(), "commit": a.commit}
    _write(out, rec)
    print(f"j7 {'within' if rec['within'] else 'beyond'}", flush=True)


# ------------------------------------------------------------------------------------------------ reproduction gate
def _gate_path(runs: Path, name: str) -> Path:
    return Path(runs) / "gate" / f"{name}.json"


def e49c_committed(name: str, path: Path = E49_READING) -> tuple[dict | None, list[str]]:
    """E49's committed per-task facts of E49c minus C0 (board, build main): ten-seed mean and SD."""
    r = json.loads(Path(path).read_text())
    if r.get("void") or not r.get("tasks") or name not in r["tasks"] or E49C not in r["tasks"][name]:
        return None, [f"{name}: E49's committed reading has no E49c line"]
    b = r["tasks"][name][E49C]["board"]
    return {"mean": float(b["mean"]), "sd": float(b["sd"]), "e49_run_commit": r.get("commit")}, []


def e49c_check(name: str, diff: np.ndarray, old: dict | None, why_old: list[str]) -> dict:
    why = list(why_old)
    worst = {"mean": None, "sd": None}
    if old is not None:
        diff = np.asarray(diff, float)
        worst = {"mean": abs(float(diff.mean()) - old["mean"]), "sd": abs(float(diff.std(ddof=1)) - old["sd"])}
        for k, v in worst.items():
            if not v <= REPRO_TOL:
                why.append(f"{name}: E49c minus C0's ten-seed {k} is {v:.4f} off E49's (tolerance {REPRO_TOL})")
    return {"worst": worst, "why": why, "pass": not why}


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    pg, why = _read_at(runs / "pregate.json", a.commit)
    if why or not pg.get("pass"):
        sys.exit("refused: the target-free gate has not passed at the run commit")
    path = _gate_path(runs, a.task)
    if path.exists():
        sys.exit(f"{path} exists")
    with e49.captured() as seen:
        cands = [harness.CANDIDATES[harness.FLOOR], e49.observed(harness.CANDIDATES[C0], seen),
                 e49.guarded(harness.CANDIDATES[E49C], seen)]
        res = harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data), build="main", log=None)
    if e49._is_test_board(seen["board"]):
        raise RuntimeError("E53 met a test board")
    scores = e47.per_seed(res)[0]
    e47_part = e49.gate_check(a.task, t.key, scores)
    old, why_old = e49c_committed(a.task)
    ov_part = e49c_check(a.task, scores[E49C] - scores[C0], old, why_old)
    rec = {"task": a.task, "key": t.key, "e47": e47_part, "e49c": {**ov_part, "e49_files_commit": E49_FILES_COMMIT},
           "e49c_refused_seeds": sorted(s for (n, s), w in seen["refused"].items() if n == E49C and w),
           "commit": res["provenance"]["commit"], "tree_clean": res["provenance"]["tree_clean"]}
    rec["pass"] = bool(e47_part["pass"] and ov_part["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
    _write(path, rec)
    print(f"gate {a.task} {'pass' if rec['pass'] else 'fail'}", flush=True)
    sys.exit(0 if rec["pass"] else 3)


def gate_problems(runs: Path, commit: str) -> list[str]:
    why = []
    pg, w = _read_at(Path(runs) / "pregate.json", commit)
    why += w
    if pg is not None and not w and not pg.get("pass"):
        why.append("the target-free gate failed")
    j, w = _read_at(Path(runs) / "j7.json", commit)
    why += w
    for name in nt.all_tasks():
        p = _gate_path(runs, name)
        if not p.exists():
            why.append(f"no gate file for {name}")
            continue
        g = json.loads(p.read_text())
        if not g.get("pass") or g.get("commit") != commit:
            why.append(f"the gate did not pass on {name} at the run commit")
    return why


# ------------------------------------------------------------------------------------------------ the run
def candidates(seen) -> list:
    return [harness.CANDIDATES[harness.FLOOR], e49.observed(harness.CANDIDATES[C0], seen),
            *(e49.guarded(harness.CANDIDATES[n], seen) for n in WRAPPED)]


def one_seed(t, s: int, data: Path) -> dict:
    """Both builds of one task-seed: run_task's results and E53's per-seed facts."""
    out: dict = {"builds": {}, "facts": {}}
    files, refused = {}, {}
    for build in BUILDS:
        with e49.captured() as seen:
            res = harness.run_task(t, candidates(seen), seeds=[s], data_dir=data, build=build, log=None)
        board = seen["board"]
        if e49._is_test_board(board):
            raise RuntimeError("E53 met a test board")
        out["builds"][build] = res
        files[build] = {n: seen["file"][(n, s)] for n in NAMES}
        refused[build] = {n: seen["refused"].get((n, s)) for n in WRAPPED}
        if build != BUILDS[0]:
            continue
        b = harness.BUILD_SEED_OFFSET + s
        recs = {n: seen["records"].get((FORM_OF.get(n), b), {}) for n in NAMES}
        truth = ev.working_copies(seen["target"], seen["ref"], s, ev.SUBSAMPLE_FRAC, harness.PROTOCOL)[0]
        out["facts"] = {"board_split": str(board.split),
                        "x": {n: seen["x"][(n, s)] for n in NAMES},
                        "n": {n: int(len(seen["stage"][(n, s)].X)) for n in NAMES},
                        "refused": refused[build], "records": recs,
                        "hands": {n: e49._hand(seen["stage"][(n, s)].coords, truth, s) for n in NAMES}}
    out["facts"]["files_identical_across_builds"] = all(files[BUILDS[0]][n] == files[b][n]
                                                        for b in BUILDS[1:] for n in NAMES)
    out["facts"]["refusals_identical_across_builds"] = all(refused[BUILDS[0]] == refused[b] for b in BUILDS[1:])
    return out


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why = gate_problems(runs, a.commit)
    if why:
        sys.exit("refused (the gates): " + "; ".join(why))
    for c in NAMES:
        if not harness.CANDIDATES[c].applies(t):
            sys.exit(f"{a.task}: {c} does not apply")
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": SCRIPT,
               **one_seed(t, s, Path(a.data))}
        _write(path, rec)
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order): E49's list."""
    why = list(why_in)
    for task in nt.all_tasks():
        rs = recs.get(task, [])
        if [r["seed"] for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
        for r in rs:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            for b in BUILDS:
                res = r["builds"].get(b)
                if res is None:
                    why.append(f"{tag}: build {b} missing")
                    continue
                p = res["provenance"]
                if p["commit"] != commit or not p["tree_clean"]:
                    why.append(f"{tag}: build {b} at another commit or on an unclean tree")
                if res["build"] != b or res["protocol"] != "bake" or [x["seed"] for x in res["rows"]] != [r["seed"]]:
                    why.append(f"{tag}: build {b} is not the bake protocol at this seed and build")
                if str(res["board"].get("split", "")).startswith("test"):
                    why.append(f"{tag}: a test board")
                row = res["rows"][0]
                if any(n not in row["scores"] or n not in row["skills"] or n not in row["groups"]
                       for n in (harness.FLOOR, *NAMES)):
                    why.append(f"{tag}: build {b} lacks an arm")
                    continue
                for n in WRAPPED:
                    d = max(abs(row["skills"][n][m] - row["skills"][C0][m]) for m in EXPR_METRICS)
                    if not d <= EXPR_TOL:
                        why.append(f"{tag}: {n}'s expression-only skills differ from C0's ({b})")
            f = r["facts"]
            if str(f.get("board_split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            if not f.get("files_identical_across_builds") or not f.get("refusals_identical_across_builds"):
                why.append(f"{tag}: a file or a refusal differs between the builds")
            jit = {n: f["records"].get(n, {}).get("n_jittered", 0) > 0 for n in NAMES}
            for n in WRAPPED:
                if f["n"][n] != f["n"][C0] or f["x"][n] != f["x"][C0]:
                    why.append(f"{tag}: {n}'s count or expression matrix differs from C0's")
                if not (jit[n] or jit[C0]) and "main" in r["builds"]:
                    row = r["builds"]["main"]["rows"][0]
                    if n in row["skills"] and C0 in row["skills"]:
                        d = abs(row["skills"][n]["scale_log_ratio"] - row["skills"][C0]["scale_log_ratio"])
                        if not d <= SCALE_TOL:
                            why.append(f"{tag}: {n}'s scale_log_ratio skill is off C0's by {d:.2e}")
    return why


def truth_means(recs, groups, a: str, b: str, what: str = "board", build: str = "main") -> dict[str, float]:
    """Each truth's mean over its tasks of the ten-seed mean difference a minus b."""
    return {truth: float(np.mean([np.mean(e49.series(recs[n], a, what, build) - e49.series(recs[n], b, what, build))
                                  for n in names]))
            for truth, names in groups.items()}


def leave_one_out(tm: dict[str, float]) -> dict[str, float]:
    """The merged mean (equal truth weights) with each truth left out in turn."""
    return {k: float(np.mean([v for j, v in tm.items() if j != k])) for k in tm}


def clauses(recs, groups, tl_groups, refused_arm: int) -> dict:
    D = {"board_e49c": e49.merged_D(recs, groups, ARM, E49C, "board"),
         "shape_c0": e49.merged_D(recs, groups, ARM, C0, f"group:{GROUP_SHAPE}"),
         "shape_e49c": e49.merged_D(recs, groups, ARM, E49C, f"group:{GROUP_SHAPE}"),
         "local_e49c": e49.merged_D(recs, groups, ARM, E49C, f"group:{GROUP_LOCAL}"),
         "S_placebo": e49.merged_D(recs, groups, ARM, PLACEBO, "S"),
         "test_like_board_e49c": e49.merged_D(recs, tl_groups, ARM, E49C, "board")}
    tm = truth_means(recs, groups, ARM, E49C)
    loo = leave_one_out(tm)
    c = {"a": bool(D["board_e49c"].mean() >= A_MERGED and sum(v > 0 for v in tm.values()) >= A_TRUTHS
                   and min(loo.values()) >= A_LEAVE_ONE_OUT),
         "b": bool(D["test_like_board_e49c"].mean() >= B_TEST_LIKE),
         "c": bool(D["shape_c0"].mean() >= C_OVER_C0 and D["shape_e49c"].mean() >= C_OVER_E49C),
         "d": bool(D["local_e49c"].mean() >= -D_LOCAL),
         "e": bool(D["S_placebo"].mean() >= E_S and (D["S_placebo"] > 0).sum() >= E_AHEAD),
         "f": bool(refused_arm <= F_MAX_REFUSED)}
    return {"clauses": c, "passes": all(c.values()), "D": {k: e49._stat(v) for k, v in D.items()},
            "truth_means_vs_e49c": tm, "truths_positive": int(sum(v > 0 for v in tm.values())),
            "leave_one_out_vs_e49c": loo, "arm_refused_task_seeds": int(refused_arm)}


def reported(recs, groups, pregate_rec: dict | None, j7_rec: dict | None) -> dict:
    per_truth = {f"{a} - {b}": truth_means(recs, groups, a, b) for a, b in
                 ((ARM, E49C), (ARM, C0), (PLACEBO, E49C), (E49C, C0))}
    per_task = {}
    for t, rs in recs.items():
        per_task[t] = {}
        for a, b in ((ARM, E49C), (ARM, C0), (ARM, PLACEBO)):
            per_task[t][f"{a} - {b}"] = {
                "board": float(np.mean(e49.series(rs, a, "board") - e49.series(rs, b, "board"))),
                **{m: float(np.mean(e49.series(rs, a, f"skill:{m}") - e49.series(rs, b, f"skill:{m}")))
                   for m in (*SHAPE_TERMS, "neighborhood_mmd")}}
        per_task[t]["arm_fit"] = {k: [r["facts"]["records"].get(ARM, {}).get(f"shape_interp_{k}") for r in rs]
                                  for k in ("lddmm_data_ratio", "jac_min", "thin_sd")}
        per_task[t]["arm_dice_proper_mirror"] = [float(np.mean([r["facts"]["hands"][ARM][k] for r in rs]))
                                                 for k in ("dice_proper", "dice_mirror")]
    return {"per_truth": per_truth, "single_truth_tasks": {t: per_task[t] for t in SINGLE_TRUTH_TASKS},
            "per_task": per_task,
            "prefix_board_vs_e49c": e49._stat(e49.merged_D(recs, groups, ARM, E49C, "board", build="prefix")),
            "refused_task_seeds": {n: int(sum(bool(r["facts"]["refused"].get(n)) for rs in recs.values() for r in rs))
                                   for n in WRAPPED},
            "gate_values": None if pregate_rec is None else pregate_rec.get("pairs"),
            "j7": j7_rec}


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str, pregate_rec: dict | None = None,
            j7_rec: dict | None = None) -> dict:
    why = check_run(recs, why_in, commit)
    if why:
        return {"void": True, "why": why, "commit": commit}
    tasks = nt.all_tasks()
    groups = e47.truths(tasks)
    tl_groups = e47.truths({n: tasks[n] for n in TEST_LIKE})
    refused_arm = sum(bool(r["facts"]["refused"].get(ARM)) for rs in recs.values() for r in rs)
    c = clauses(recs, groups, tl_groups, refused_arm)
    outcome = ("all of (a) to (f) hold: the candidate file is fixed in the decision log before any score of it and "
               "goes to the captain as a priced question beside E49c" if c["passes"] else
               "not all of (a) to (f) hold: no file")
    return {"void": False, "why": [], "commit": commit, "tries": TRY_IDS, "truths": groups,
            "test_like_truths": tl_groups, **c, "outcome": outcome,
            "reported": reported(recs, groups, pregate_rec, j7_rec)}


def text(res: dict) -> str:
    if res["void"]:
        return f"Row E53 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n"
    D = res["D"]
    L = [f"Row E53 (run commit {res['commit']}): E53a {ARM} (the arm), E53b {PLACEBO} (the placebo); the target-free "
         f"gate and the reproduction gate passed (C0 against E47, E49c against E49, {REPRO_TOL} board points).",
         "Pass line, build main, merged over five truths (mean, SE, seeds ahead of 10):",
         f"  (a) {'holds' if res['clauses']['a'] else 'fails'}: board vs E49c {D['board_e49c']['mean']:+.3f} (SE "
         f"{D['board_e49c']['se']:.3f}); truths positive {res['truths_positive']} of 5; worst leave-one-out "
         f"{min(res['leave_one_out_vs_e49c'].values()):+.3f}",
         f"  (b) {'holds' if res['clauses']['b'] else 'fails'}: test-like board vs E49c "
         f"{D['test_like_board_e49c']['mean']:+.3f}",
         f"  (c) {'holds' if res['clauses']['c'] else 'fails'}: gated shape_growth vs C0 {D['shape_c0']['mean']:+.4f}, "
         f"vs E49c {D['shape_e49c']['mean']:+.4f}",
         f"  (d) {'holds' if res['clauses']['d'] else 'fails'}: local_organisation vs E49c "
         f"{D['local_e49c']['mean']:+.4f}",
         f"  (e) {'holds' if res['clauses']['e'] else 'fails'}: ungated S vs the placebo {D['S_placebo']['mean']:+.4f} "
         f"({D['S_placebo']['ahead']} seeds ahead)",
         f"  (f) {'holds' if res['clauses']['f'] else 'fails'}: the arm refused on {res['arm_refused_task_seeds']} of "
         f"90 task-seeds",
         f"Outcome: {res['outcome']}.",
         "Reported, deciding nothing:",
         "  per truth, arm minus E49c board: " + ", ".join(f"{k} {v:+.3f}" for k, v in
                                                            res["reported"]["per_truth"][f"{ARM} - {E49C}"].items())]
    for t, v in res["reported"]["single_truth_tasks"].items():
        x = v[f"{ARM} - {E49C}"]
        L.append(f"  {t}: arm minus E49c board {x['board']:+.3f}, d2_shape {x['d2_shape']:+.4f}, occupancy_dice "
                 f"{x['occupancy_dice']:+.4f}")
    L.append(f"  prefix build, board vs E49c {res['reported']['prefix_board_vs_e49c']['mean']:+.3f}; refusals "
             + ", ".join(f"{k} {v}" for k, v in res["reported"]["refused_task_seeds"].items()))
    j = res["reported"].get("j7") or {}
    if j:
        L.append(f"  J7 dry run on E48's committed rescore: worst {j['worst']['seed']:.4f} per seed, "
                 f"{j['worst']['mean']:.4f} per mean, tolerance {j['tolerance']} "
                 f"({'within' if j['within'] else 'beyond'})")
    L.append("Both tries are diagnostics (J5). Nothing here shows E7.75's shape, the live gate, build or metric "
             "seed.")
    return "\n".join(L) + "\n"


def provenance() -> dict:
    s = repo_state(REPO)
    return {"head": s.commit, "clean": s.unclean == (), "script_sha256": e49._sha_file(Path(__file__)),
            "shape_interp_sha256": e49._sha_file(Path(si.__file__)),
            "lddmm_sha256": e49._sha_file(Path(lddmm.__file__))}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    if not has_row():
        sys.exit("E53 has no row in experiments/harness/TRIES.md")
    runs = Path(a.runs)
    why = [*gate_problems(runs, a.commit), *commit_problems(a.commit)]
    pg = _read_at(runs / "pregate.json", a.commit)[0]
    jr = _read_at(runs / "j7.json", a.commit)[0]
    res = {**reading(e49.load_runs(runs), why, a.commit, pg, jr), "reading_code": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(text(res))
    print(f"wrote {files[1]}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("describe")
    d.add_argument("--pair", required=True, help="'test' or a short name of nine_tasks.all_tasks()")
    for name in ("pregate", "j7", "read"):
        sub.add_parser(name)
    for name in ("gate", "run"):
        p = sub.add_parser(name)
        p.add_argument("--task", required=True, help="a short name of nine_tasks.all_tasks()")
    for p in sub.choices.values():
        p.add_argument("--commit", required=True, help="the run commit (full hash)")
        p.add_argument("--runs", default=str(RUNS))
        p.add_argument("--data", default=str(harness.DATA_DIR))
    sub.choices["read"].add_argument("--out", default=str(RESULTS))
    a = ap.parse_args(argv)
    {"describe": cmd_describe, "pregate": cmd_pregate, "j7": cmd_j7, "gate": cmd_gate, "run": cmd_run,
     "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
