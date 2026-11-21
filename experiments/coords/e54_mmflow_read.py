#!/usr/bin/env python
"""Row E54 of experiments/harness/TRIES.md (the moonshot "multi-stage flow, then retrieval of real cells", specified by
a score-blind judge; the specification's own labels E52a-E52c are E54a-E54c here; the row is placed by the main session
before any run). The kill test. No number here comes from a test score; the script meets no test board.

Arms. E54b fgw_mmflow_snap (the arm) and E54c fgw_mmflow_anti (its placebo), harness candidates made around one task's
fitted states (vec.harness.mmflow_candidate): fgw_shape0x's positions, count, bracket quotas and coordinates, every
slot's cell replaced by the real cell of its own bracket retrieved one to one (vec.t2_mmflow.retrieve, rectangular,
against the bracket's whole pool) for the state its drawn cell reaches under the flow: integrated to the target day
(snap) or away from it by the same span (anti). C0 = fgw_shape0x and the floor are rescored and not counted.

Fits (vec.t2_mmflow, constants in its header). For each practice task the embedding and the flow are fitted without
the target stage and without every released stage strictly inside the task's bracket (EXCLUSIONS below; seven distinct
fits), and the gate adds the leave-one-out fit of E8.25_late: eight fits (FOLDS). The training stages are the seven
released stages of the embryo series (E6.75, E7.25, E7.5, E8.0 and the heart-setting files E8.25_late, E8.5, E8.75)
less a fit's exclusions; the genes the validation board's embryo panel (every one of the seven files holds it: run_task
refuses an input that lacks a panel gene, and each of the seven is an input of a practice task).

Steps (each refuses unless row E54 is in TRIES.md, the tree is clean and HEAD is the run commit):
  fit --fold NAME --commit SHA    target-free (it loads only the fold's training stages; no try): the embedding, the
                                  flow, and the integrated states every task and gate stage of the fold needs, under
                                  runs/e54/fits/<fold>/. A fold whose done.json is of this commit is not refitted.
  repro --task NAME --commit SHA  the reproduction gate (not counted): the floor and C0 under build main, ten seeds,
                                  against E47's committed scores (run at 9967cf1, committed in d5ee9b5) within 0.01
                                  board points per seed and per ten-seed mean (E49's check, J7's tolerance); exit 3 on
                                  a miss. A file of this commit is read again, not re-run.
  gate --commit SHA               E54a (counted, J5 diagnostic): leave-one-stage-out on E7.25, E7.5, E8.0, E8.25_late
                                  and E8.5, each in its fit's own thinned latent: energy distance to 5 draws of 5,000
                                  held cells of the flow-predicted states (both neighbours carried to the day) against
                                  the calendar-fraction mixture of the same neighbour cells as measured. Passes when the
                                  flow is lower on at least 4 of 5 stages and on E7.5. Writes the gate's reading
                                  (gate.json, gate.txt under the results folder, never overwritten); exit 3 if it fails,
                                  and then E54b and E54c do not run and are not counted.
  run --task NAME --commit SHA    E54b and E54c with C0 and the floor, build main, bake protocol, seeds 0-9, build seed
                                  1,000,000 + seed, seed by seed (resumable); refused unless the gate passed, every
                                  repro file passes and the task's fit is done, all at this commit
  read --commit SHA               the kill test's reading (reading.json, reading.txt), never overwritten
The script prints only "fit <fold> done", "repro <task> pass|fail", "gate pass|fail", "<task> seed <s> done" and
"wrote <path>". It prints no count.

Reading, fixed before the run (applied by `reading`). D_s(X) for a difference X at seed s is the mean over the five
truths of each truth's task-mean difference at s (E8.0: half, late; E8.25_late: new1, x1; E8.5: new2, x2, x3; E7.5:
e75; E7.25: e725); merged = mean over s. Test-like: half, late, new1, new2 (three truths). Under build main the kill
test passes only if all hold:
  (a) merged board gain of the arm over C0 at least +1.0; the per-truth mean gains positive on at least 4 of 5 truths;
      and with any one truth left out the merged gain over the other four at least +0.5;
  (b) test-like merged board gain over C0 at least +1.0;
  (c) arm minus anti merged board at least +0.7, with D_s above 0 in at least 8 of 10 seeds;
  (d) merged gated cell_state group skill gain over C0 above 0;
  (e) the merged board gain without e725 (the other four truths) at least +0.7.
Void (read as nothing): either bracket's median library size (vec.t2_build.library_sizes over the drawn rows) of the
arm or the anti more than 5% off C0's at any task-seed; a shape term (d2_shape, occupancy_dice, scale_log_ratio skill)
off C0's by more than 1e-9 where neither file was jittered; a count differing from C0's; a missing task, seed or arm
(a refused arm included); results from another commit or an unclean tree; a repro file missing or failing; the gate
reading missing, from another commit, or not passed; a test board.
Reported, deciding nothing: per-stage gate distances; each bracket's Kish pool use; the type shares drawn against C0's;
e75 and e725 apart; every truth's gain.
Counts: E54a always once the gate has run; E54b and E54c only if the gate passed. All are J5 diagnostics.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))
sys.path.insert(0, str(REPO / "experiments" / "harness"))
sys.path.insert(0, str(REPO / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
import e47_floor_exact as e47  # noqa: E402
import e49_shape_interp_read as e49  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import evaluate as ev  # noqa: E402
from vec import harness  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec import t2_mmflow as mf  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW_IDS = ("E54", "E54a")  # the row's first word: "E54" as the brief labels it, or "E54a to E54c"
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
TRY_IDS = {"gate": "E54a", "fgw_mmflow_snap": "E54b", "fgw_mmflow_anti": "E54c"}
C0 = "fgw_shape0x"
ARM, ANTI = "fgw_mmflow_snap", "fgw_mmflow_anti"
ARMS = (ARM, ANTI)
STAGE_FILES = ("T2/E6.75.h5ad", "T2/E7.25.h5ad", "T2/E7.5.h5ad", "T2/E8.0.h5ad", "T2/E8.25_late.h5ad",
               "T2/E8.5.h5ad", "T2/E8.75.h5ad")
DAYS = {f: float(harness.STAGE_DAYS[f]) for f in STAGE_FILES}
PANEL_BOARD = "T2:embryo:val_interp"
# the specification's exclusion sets (its stage 8.25 is the file E8.25_late, day 8.3 in vec.harness.STAGE_DAYS)
EXCLUSIONS = {"half": ("T2/E8.0.h5ad", "T2/E8.25_late.h5ad"), "late": ("T2/E8.0.h5ad",),
              "new1": ("T2/E8.0.h5ad", "T2/E8.25_late.h5ad"),
              "new2": ("T2/E8.0.h5ad", "T2/E8.25_late.h5ad", "T2/E8.5.h5ad"),
              "x1": ("T2/E8.25_late.h5ad", "T2/E8.5.h5ad"), "x2": ("T2/E8.25_late.h5ad", "T2/E8.5.h5ad"),
              "x3": ("T2/E8.5.h5ad",), "e75": ("T2/E7.5.h5ad",), "e725": ("T2/E7.25.h5ad",)}
GATE_STAGES = ("T2/E7.25.h5ad", "T2/E7.5.h5ad", "T2/E8.0.h5ad", "T2/E8.25_late.h5ad", "T2/E8.5.h5ad")
GATE_MUST = "T2/E7.5.h5ad"
GATE_NEED = 4
TEST_LIKE = ("half", "late", "new1", "new2")
DROP_E725 = "e725"
BAR_A, BAR_A_LOO, BAR_B, BAR_C, BAR_E = 1.0, 0.5, 1.0, 0.7, 0.7
TRUTHS_POSITIVE = 4
AHEAD_C = 8
LIB_TOL = 0.05
SHAPE_TERMS = ("d2_shape", "occupancy_dice", "scale_log_ratio")
SHAPE_TOL = 1e-9
GROUP_CELL_STATE = "cell_state"
BASES = ("c902b90",)  # main when the branch was cut
E47_DIR = e49.E47_DIR
E47_FILES_COMMIT = e49.E47_FILES_COMMIT
SCRIPT = "experiments/coords/e54_mmflow_read.py"
RUNS = REPO / "runs" / "e54"
RESULTS = REPO / "experiments" / "coords" / "results" / "e54_mmflow"


def _stem(f: str) -> str:
    return Path(f).stem.removeprefix("E").replace(".", "")


def fold_name(excluded) -> str:
    return "x" + "_".join(_stem(f) for f in sorted(excluded, key=DAYS.get))


def task_folds() -> dict[str, str]:
    return {name: fold_name(ex) for name, ex in EXCLUSIONS.items()}


def gate_folds() -> dict[str, str]:
    return {f: fold_name((f,)) for f in GATE_STAGES}


def folds() -> dict[str, tuple[str, ...]]:
    """Fold name -> the stages it leaves out: the tasks' seven and the gate's one more."""
    out = {fold_name(ex): tuple(sorted(ex, key=DAYS.get)) for ex in EXCLUSIONS.values()}
    out.update({fold_name((f,)): (f,) for f in GATE_STAGES})
    return out


def exclusion_problems() -> list[str]:
    """Each task's exclusions as the rule makes them (its target and every released stage strictly inside its bracket)
    against the specification's table."""
    why = []
    for name, t in nt.all_tasks().items():
        rule = mf.exclusion(DAYS[t.inputs["before"]], DAYS[t.inputs["after"]], DAYS[t.target], list(DAYS.values()))
        if rule != frozenset(DAYS[f] for f in EXCLUSIONS[name]):
            why.append(f"{name}: the rule's exclusions differ from the specification's")
        if t.inputs["before"] in EXCLUSIONS[name] or t.inputs["after"] in EXCLUSIONS[name]:
            why.append(f"{name}: a bracket is excluded")
    return why


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def has_row(path: Path | None = None) -> bool:
    rows = n12.tries_rows(TRIES if path is None else path)
    return any(r in rows for r in ROW_IDS)


def commit_problems(commit: str) -> list[str]:
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    for base in BASES:
        if _git("merge-base", "--is-ancestor", base, "HEAD").returncode != 0:
            why.append(f"HEAD does not contain {base}")
    if _git("diff", "--quiet", E47_FILES_COMMIT, "HEAD", "--", str(E47_DIR.relative_to(REPO))).returncode != 0:
        why.append(f"E47's committed results differ from {E47_FILES_COMMIT}'s")
    return why


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (ev, harness, tb, mf):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    if not has_row():
        sys.exit("E54 has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E54 runs on a committed tree")
    why = commit_problems(commit) + exclusion_problems()
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E54 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _write_new(path: Path, text: str) -> None:
    """Write a file that must not exist yet (a reading is never overwritten)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "x", encoding="utf-8") as f:
        f.write(text)


def _write_replace(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


# ------------------------------------------------------------------------------------------------ fit (target-free)
def _panel() -> list[str]:
    from vec.boards import load_boards

    board = load_boards()[PANEL_BOARD]
    if _is_test_board(board):
        sys.exit("the panel comes from the validation board only")
    return list(board.genes)


def latent_space(stages: dict[str, object]) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """The fit's embedding from its training stages (file -> loaded AnnData on the panel): thinned features, the PCA,
    and every training cell's latent."""
    feats = {DAYS[f]: mf.features(mf.thin_counts(a.layers["counts"], DAYS[f])) for f, a in stages.items()}
    centre, loadings = mf.fit_embedding(feats)
    return centre, loadings, {f: mf.embed(feats[DAYS[f]], centre, loadings) for f in stages}


def project(a, day: float, centre: np.ndarray, loadings: np.ndarray) -> np.ndarray:
    """A stage the fit left out, thinned on its own stream and projected on the fit's PCA (the gate's held stage)."""
    return mf.embed(mf.features(mf.thin_counts(a.layers["counts"], day)), centre, loadings)


def fold_outputs(fold: str, stages: dict[str, object], lat: dict[str, np.ndarray], net, device: str) -> dict:
    """The integrated states of every task and gate stage of the fold (file name -> arrays)."""
    out = {}
    for name, fname in task_folds().items():
        if fname != fold:
            continue
        t = nt.all_tasks()[name]
        b, a, tgt = DAYS[t.inputs["before"]], DAYS[t.inputs["after"]], DAYS[t.target]
        Lb, La = lat[t.inputs["before"]], lat[t.inputs["after"]]
        out[f"states_{name}.npz"] = {
            "before_obs": np.asarray(stages[t.inputs["before"]].obs_names, dtype=str),
            "after_obs": np.asarray(stages[t.inputs["after"]].obs_names, dtype=str),
            "before_latent": Lb, "after_latent": La,
            "before_snap": mf.integrate(net, Lb, b, tgt, device), "after_snap": mf.integrate(net, La, a, tgt, device),
            "before_anti": mf.integrate(net, Lb, b, b - (tgt - b), device),
            "after_anti": mf.integrate(net, La, a, a + (a - tgt), device),
            "days": np.array([b, a, tgt], dtype=np.float64)}
    for f, fname in gate_folds().items():
        if fname != fold:
            continue
        lo, hi = mf.neighbours(DAYS[f], list(DAYS.values()))
        flo, fhi = (next(k for k, d in DAYS.items() if d == x) for x in (lo, hi))
        out[f"gate_{_stem(f)}.npz"] = {
            "lo_latent": lat[flo], "hi_latent": lat[fhi],
            "lo_fwd": mf.integrate(net, lat[flo], lo, DAYS[f], device),
            "hi_bwd": mf.integrate(net, lat[fhi], hi, DAYS[f], device),
            "days": np.array([lo, hi, DAYS[f]], dtype=np.float64)}
    return out


def cmd_fit(a) -> None:
    _refuse_unless_ready(a.commit)
    import torch

    fold = a.fold
    if fold not in folds():
        sys.exit(f"unknown fold {fold}; known: {sorted(folds())}")
    folder = Path(a.runs) / "fits" / fold
    done = folder / "done.json"
    if done.exists():
        rec = json.loads(done.read_text())
        if rec.get("commit") != a.commit:
            sys.exit(f"{done} exists from another commit")
        print(f"fit {fold} done", flush=True)
        return
    if a.device.startswith("cuda") and not torch.cuda.is_available():
        sys.exit("no GPU visible: the fits run on a GPU (CUDA_VISIBLE_DEVICES)")
    from vec.prep import load_stage

    genes = _panel()
    included = [f for f in STAGE_FILES if f not in folds()[fold]]
    stages = {f: load_stage(Path(a.data) / f, genes=genes, drop_duplicates=True) for f in included}
    centre, loadings, lat = latent_space(stages)
    net, losses = mf.fit_flow([lat[f] for f in included], [DAYS[f] for f in included], device=a.device)
    folder.mkdir(parents=True, exist_ok=True)
    torch.save(net.state_dict(), folder / "model.pt")
    np.savez(folder / "embed.npz", centre=centre, loadings=loadings, genes=np.asarray(genes, dtype=str),
             included=np.asarray(included, dtype=str))
    written = {}
    for name, arrays in fold_outputs(fold, stages, lat, net, a.device).items():
        np.savez(folder / name, **arrays)
        written[name] = _sha_file(folder / name)
    rec = {"commit": a.commit, "fold": fold, "excluded": list(folds()[fold]), "included": included,
           "rule": mf.RULE, "device": a.device, "torch": torch.__version__,
           "loss_first_100": float(losses[:100].mean()), "loss_last_100": float(losses[-100:].mean()),
           "files_sha256": {"model.pt": _sha_file(folder / "model.pt"), "embed.npz": _sha_file(folder / "embed.npz"),
                            **written}}
    _write_replace(done, json.dumps(rec, indent=1))
    print(f"fit {fold} done", flush=True)


def fit_problems(runs: Path, commit: str, wanted) -> list[str]:
    why = []
    for fold in sorted(set(wanted)):
        p = Path(runs) / "fits" / fold / "done.json"
        if not p.exists():
            why.append(f"fold {fold} is not fitted")
            continue
        rec = json.loads(p.read_text())
        if rec.get("commit") != commit:
            why.append(f"fold {fold} was fitted at another commit")
            continue
        for name, sha in rec.get("files_sha256", {}).items():
            f = p.parent / name
            if not f.exists() or _sha_file(f) != sha:
                why.append(f"fold {fold}: {name} differs from its fit")
    return why


def load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


# ------------------------------------------------------------------------------------------------ repro (not counted)
def _repro_path(runs: Path, name: str) -> Path:
    return Path(runs) / "repro" / f"{name}.json"


def cmd_repro(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _repro_path(Path(a.runs), a.task)
    if path.exists():
        rec = json.loads(path.read_text())
        if rec.get("commit") != a.commit:
            sys.exit(f"{path} exists from another commit")
    else:
        res = harness.run_task(t, [harness.CANDIDATES[harness.FLOOR], harness.CANDIDATES[C0]], seeds=SEEDS,
                               data_dir=Path(a.data), build="main", log=None)
        if str(res["board"].get("split", "")).startswith("test"):
            sys.exit("E54 met a test board")
        rec = {**e49.gate_check(a.task, t.key, e47.per_seed(res)[0]), "commit": res["provenance"]["commit"],
               "tree_clean": res["provenance"]["tree_clean"]}
        rec["pass"] = bool(rec["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
        _write_replace(path, json.dumps(rec, indent=1, default=float))
    print(f"repro {a.task} {'pass' if rec['pass'] else 'fail'}", flush=True)
    sys.exit(0 if rec["pass"] else 3)


def repro_problems(runs: Path, commit: str) -> list[str]:
    why = []
    for name in nt.all_tasks():
        p = _repro_path(runs, name)
        if not p.exists():
            why.append(f"no repro file for {name}")
            continue
        g = json.loads(p.read_text())
        if not g.get("pass") or g.get("commit") != commit:
            why.append(f"the reproduction gate did not pass on {name} at the run commit")
    return why


# ------------------------------------------------------------------------------------------------ gate (E54a)
def gate_reading(stages: dict[str, dict], commit: str) -> dict:
    passed = mf.gate_pass(stages, GATE_MUST, GATE_NEED)
    return {"try": TRY_IDS["gate"], "commit": commit, "stages": stages, "must": GATE_MUST, "need": GATE_NEED,
            "beats": int(sum(bool(v["beats"]) for v in stages.values())), "pass": passed,
            "outcome": ("the gate passes: E54b and E54c run and are counted" if passed else
                        "the gate fails: E54b and E54c do not run and are not counted; E54a stays counted"),
            "counted": [TRY_IDS["gate"], TRY_IDS[ARM], TRY_IDS[ANTI]] if passed else [TRY_IDS["gate"]]}


def gate_text(g: dict) -> str:
    L = [f"Row E54, {g['try']} (run commit {g['commit']}): leave-one-stage-out, energy distance in each fit's own "
         "thinned latent, mean of 5 draws of 5,000 held cells (flow-predicted against the calendar-fraction mixture):"]
    for f, v in g["stages"].items():
        L.append(f"  {f}: flow {v['flow_mean']:.5f}, mixture {v['mixture_mean']:.5f}; the flow "
                 f"{'beats' if v['beats'] else 'does not beat'} the mixture")
    L.append(f"The flow beats the mixture on {g['beats']} of {len(g['stages'])} stages (needed: {g['need']} and "
             f"{g['must']}). Outcome: {g['outcome']}. Counted: {', '.join(g['counted'])} (J5 diagnostics).")
    return "\n".join(L) + "\n"


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    out = Path(a.out)
    files = [out / "gate.json", out / "gate.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    runs = Path(a.runs)
    why = fit_problems(runs, a.commit, gate_folds().values()) + repro_problems(runs, a.commit)
    if why:
        sys.exit("refused (the fits and the reproduction gate): " + "; ".join(why))
    from vec.prep import load_stage

    genes = _panel()
    stages = {}
    for f, fold in gate_folds().items():
        emb = load_npz(runs / "fits" / fold / "embed.npz")
        if list(emb["genes"]) != genes:
            sys.exit(f"fold {fold} was fitted on another panel")
        st = load_npz(runs / "fits" / fold / f"gate_{_stem(f)}.npz")
        lo, hi, day = (float(x) for x in st["days"])
        held = project(load_stage(Path(a.data) / f, genes=genes, drop_duplicates=True), day, emb["centre"],
                       emb["loadings"])
        stages[f] = {"neighbours": [lo, hi], "fraction": (day - lo) / (hi - lo),
                     **mf.gate_stage(st["lo_latent"], st["lo_fwd"], st["hi_latent"], st["hi_bwd"], held,
                                     (day - lo) / (hi - lo))}
    g = {**gate_reading(stages, a.commit), "reading_code": provenance()}
    _write_new(files[0], json.dumps(g, indent=1, default=float))
    _write_new(files[1], gate_text(g))
    print(f"gate {'pass' if g['pass'] else 'fail'}", flush=True)
    sys.exit(0 if g["pass"] else 3)


def load_gate(out: Path, commit: str) -> tuple[dict | None, list[str]]:
    p = Path(out) / "gate.json"
    if not p.exists():
        return None, ["no gate reading (E54a)"]
    g = json.loads(p.read_text())
    why = [] if g.get("commit") == commit else ["the gate reading is from another commit"]
    if not g.get("pass"):
        why.append("the gate did not pass: E54b and E54c are not run")
    return g, why


# ------------------------------------------------------------------------------------------------ run (E54b, E54c)
def file_digest(stage) -> str:
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    C = np.ascontiguousarray(stage.coords, dtype=np.float32)
    return hashlib.sha256(np.asarray(X.shape, np.int64).tobytes() + X.tobytes() + C.tobytes()).hexdigest()


def coords_digest(stage) -> str:
    return hashlib.sha256(np.ascontiguousarray(stage.coords, dtype=np.float32).tobytes()).hexdigest()


@contextlib.contextmanager
def captured():
    seen: dict = {"stage": {}, "refused": {}}
    orig = harness.score_seeds

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(task=t, board=board)
        return orig(t, board, target, ref, inputs, candidates, *a, **kw)

    harness.score_seeds = score_seeds
    try:
        yield seen
    finally:
        harness.score_seeds = orig


def observed(base, seen):
    def build(ctx):
        pred = base.build(ctx)
        seen["stage"][(base.name, ctx.seed - harness.BUILD_SEED_OFFSET)] = pred
        return pred

    return dataclasses.replace(base, build=build)


def guarded(base, seen):
    """An E54 arm: its file, or on a refusal (vec.t2_mmflow.MmflowRefused) C0's file of the seed, recorded (and the
    reading voids: a missing arm)."""
    def build(ctx):
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        if _is_test_board(seen["board"]):
            raise RuntimeError("E54 met a test board")
        try:
            pred, why = base.build(ctx), None
        except mf.MmflowRefused as e:
            pred, why = seen["stage"][(C0, s)], str(e)
        seen["refused"][(base.name, s)] = why
        seen["stage"][(base.name, s)] = pred
        return pred

    return dataclasses.replace(base, build=build)


def candidates(states: dict, seen, facts: dict) -> list:
    return [harness.CANDIDATES[harness.FLOOR], observed(harness.CANDIDATES[C0], seen),
            *(guarded(harness.mmflow_candidate(n, states, facts[n]), seen) for n in ARMS)]


def one_seed(t, s: int, data: Path, states: dict) -> dict:
    facts = {n: {} for n in ARMS}
    with captured() as seen:
        res = harness.run_task(t, candidates(states, seen, facts), seeds=[s], data_dir=data, build="main", log=None)
    if _is_test_board(seen["board"]):
        raise RuntimeError("E54 met a test board")
    b = harness.BUILD_SEED_OFFSET + s
    c0 = seen["stage"][(C0, s)]
    return {"builds": {"main": res},
            "facts": {"board_split": str(seen["board"].split),
                      "refused": {n: seen["refused"].get((n, s)) for n in ARMS},
                      "arm": {n: facts[n].get(b) for n in ARMS},
                      "coords_equal_c0": {n: coords_digest(seen["stage"][(n, s)]) == coords_digest(c0) for n in ARMS},
                      "file_differs_from_c0": {n: file_digest(seen["stage"][(n, s)]) != file_digest(c0)
                                               for n in ARMS}}}


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    _, why = load_gate(Path(a.out), a.commit)
    why += repro_problems(runs, a.commit) + fit_problems(runs, a.commit, [task_folds()[a.task]])
    if why:
        sys.exit("refused (the gate, the reproduction gate and the fit): " + "; ".join(why))
    if not harness.CANDIDATES[C0].applies(t):
        sys.exit(f"{a.task}: {C0} does not apply")
    states = load_npz(runs / "fits" / task_folds()[a.task] / f"states_{a.task}.npz")
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        res = one_seed(t, s, Path(a.data), states)
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": SCRIPT,
               "fold": task_folds()[a.task], **res}
        _write_replace(path, json.dumps(rec, indent=1, default=float))
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def _row(r: dict) -> dict:
    return r["builds"]["main"]["rows"][0]


def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order)."""
    why = list(why_in)
    names = (harness.FLOOR, C0, *ARMS)
    for task in nt.all_tasks():
        rs = recs.get(task, [])
        if [r["seed"] for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
        for r in rs:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            res = r["builds"].get("main")
            if res is None:
                why.append(f"{tag}: build main missing")
                continue
            p = res["provenance"]
            if p["commit"] != commit or not p["tree_clean"]:
                why.append(f"{tag}: run at another commit or on an unclean tree")
            if res["build"] != "main" or res["protocol"] != "bake" or [x["seed"] for x in res["rows"]] != [r["seed"]]:
                why.append(f"{tag}: not the bake protocol at this seed under build main")
            f = r["facts"]
            if str(res["board"].get("split", "")).startswith("test") or str(f.get("board_split", "")).startswith(
                    "test"):
                why.append(f"{tag}: a test board")
            row = res["rows"][0]
            if any(n not in row["scores"] or n not in row["skills"] for n in names):
                why.append(f"{tag}: an arm is missing")
                continue
            for n in ARMS:
                if f["refused"].get(n):
                    why.append(f"{tag}: {n} refused ({f['refused'][n]}): a missing arm")
                    continue
                if row["n_cells"][n] != row["n_cells"][C0]:
                    why.append(f"{tag}: {n}'s count differs from C0's")
                rec = f["arm"].get(n) or {}
                if not rec.get("n_jittered_any", False):
                    d = max(abs(row["skills"][n][m] - row["skills"][C0][m]) for m in SHAPE_TERMS)
                    if not d <= SHAPE_TOL:
                        why.append(f"{tag}: {n}'s shape terms are off C0's by {d:.2e}")
                for br in ("before", "after"):
                    ratio = (rec.get(br) or {}).get("library_ratio")
                    if ratio is None or not abs(ratio - 1.0) <= LIB_TOL:
                        why.append(f"{tag}: {n}'s {br} bracket median library is "
                                   + ("unreadable" if ratio is None else f"{ratio:.3f} of C0's"))
    return why


def series(recs: list[dict], name: str, what: str) -> np.ndarray:
    rows = [_row(r) for r in recs]
    if what == "board":
        return np.array([x["scores"][name] for x in rows], float)
    if what.startswith("group:"):
        return np.array([x["groups"][name][what[6:]] for x in rows], float)
    if what.startswith("skill:"):
        return np.array([x["skills"][name][what[6:]] for x in rows], float)
    raise ValueError(what)


def merged_D(recs, groups, a: str, b: str, what: str = "board") -> np.ndarray:
    return e47.merged_seeds({n: series(recs[n], a, what) - series(recs[n], b, what)
                             for names in groups.values() for n in names}, groups)


def per_truth(recs, groups, a: str, b: str, what: str = "board") -> dict[str, float]:
    return {truth: float(np.mean([np.mean(series(recs[n], a, what) - series(recs[n], b, what)) for n in names]))
            for truth, names in groups.items()}


def _stat(D: np.ndarray) -> dict:
    return {"mean": float(D.mean()), "se": float(D.std(ddof=1) / np.sqrt(len(D))) if len(D) > 1 else None,
            "ahead": int((D > 0).sum()), "behind": int((D < 0).sum())}


def pass_line(recs: dict[str, list[dict]]) -> dict:
    tasks = nt.all_tasks()
    groups = e47.truths(tasks)
    tl = e47.truths({n: tasks[n] for n in TEST_LIKE})
    no725 = e47.truths({n: t for n, t in tasks.items() if n != DROP_E725})
    board = merged_D(recs, groups, ARM, C0)
    truths = per_truth(recs, groups, ARM, C0)
    loo = {left: float(merged_D(recs, {k: v for k, v in groups.items() if k != left}, ARM, C0).mean())
           for left in groups}
    test_like = merged_D(recs, tl, ARM, C0)
    vs_anti = merged_D(recs, groups, ARM, ANTI)
    cell_state = merged_D(recs, groups, ARM, C0, f"group:{GROUP_CELL_STATE}")
    without = merged_D(recs, no725, ARM, C0)
    clauses = {
        "a": bool(board.mean() >= BAR_A and sum(v > 0 for v in truths.values()) >= TRUTHS_POSITIVE
                  and all(v >= BAR_A_LOO for v in loo.values())),
        "b": bool(test_like.mean() >= BAR_B),
        "c": bool(vs_anti.mean() >= BAR_C and (vs_anti > 0).sum() >= AHEAD_C),
        "d": bool(cell_state.mean() > 0),
        "e": bool(without.mean() >= BAR_E)}
    return {"clauses": clauses, "passes": all(clauses.values()),
            "merged_board": _stat(board), "per_truth": truths, "leave_one_truth_out": loo,
            "test_like_board": _stat(test_like), "arm_minus_anti": _stat(vs_anti),
            "cell_state": _stat(cell_state), "without_e725": _stat(without),
            "truths": groups, "test_like_truths": tl}


def _mean_or_none(values) -> float | None:
    vals = [v for v in values if v is not None]
    return float(np.mean(vals)) if vals else None


def reported(recs: dict[str, list[dict]]) -> dict:
    """Deciding nothing: per task, the arm's and the anti's board gain over C0, each bracket's Kish pool use and
    library ratio (means over seeds), and the type shares drawn against C0's (mean over seeds)."""
    out = {}
    for task, rs in recs.items():
        per = {}
        for n in ARMS:
            d = series(rs, n, "board") - series(rs, C0, "board")
            br = {}
            for b in ("before", "after"):
                facts = [(r["facts"]["arm"].get(n) or {}).get(b) or {} for r in rs]
                shares_arm, shares_c0 = {}, {}
                for f in facts:
                    for k, v in f.get("shares_arm", {}).items():
                        shares_arm[k] = shares_arm.get(k, 0.0) + v / len(facts)
                    for k, v in f.get("shares_c0", {}).items():
                        shares_c0[k] = shares_c0.get(k, 0.0) + v / len(facts)
                br[b] = {key: _mean_or_none([f.get(key) for f in facts])
                         for key in ("kish_frac_c0", "kish_frac_arm", "library_ratio")}
                br[b]["shares_arm"], br[b]["shares_c0"] = shares_arm, shares_c0
            per[n] = {"board_gain": float(d.mean()), "ahead": int((d > 0).sum()), "brackets": br}
        out[task] = per
    return out


def reading(recs: dict[str, list[dict]], gate: dict | None, why_in: list[str], commit: str) -> dict:
    why = check_run(recs, why_in, commit)
    if gate is None or not gate.get("pass"):
        why.append("no passed gate reading: E54b and E54c are not read")
    if why:
        return {"void": True, "why": why, "commit": commit, "gate": gate}
    p = pass_line(recs)
    rep = reported(recs)
    outcome = ("the kill test passes: (a) to (e) all hold" if p["passes"] else
               "the kill test fails: " + ", ".join(f"({k}) fails" for k, v in p["clauses"].items() if not v))
    return {"void": False, "why": [], "commit": commit, "gate": gate, "pass_line": p, "passes": p["passes"],
            "outcome": outcome, "counted": [TRY_IDS["gate"], TRY_IDS[ARM], TRY_IDS[ANTI]], "reported": rep,
            "e75_e725_apart": {n: {k: rep[n][a]["board_gain"] for a, k in ((ARM, "arm"), (ANTI, "anti"))}
                               for n in ("e75", "e725")}}


def _fmt(x, f="+.3f") -> str:
    return "n/a" if x is None else format(x, f)


def text(res: dict) -> str:
    if res["void"]:
        return f"Row E54 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n"
    p = res["pass_line"]
    m = p["merged_board"]
    L = [f"Row E54 (run commit {res['commit']}): gate E54a passed ({res['gate']['beats']} of "
         f"{len(res['gate']['stages'])} stages and {res['gate']['must']}); arm {ARM} (E54b), placebo {ANTI} (E54c), "
         f"against C0 = {C0}, build main, merged over five truths.",
         f"(a) merged board gain {m['mean']:+.3f} (SE {_fmt(m['se'], '.3f')}, {m['ahead']} of 10 ahead); per truth "
         + ", ".join(f"{k} {v:+.3f}" for k, v in p["per_truth"].items()) + "; leaving one truth out "
         + ", ".join(f"{k} {v:+.3f}" for k, v in p["leave_one_truth_out"].items())
         + f": {'holds' if p['clauses']['a'] else 'fails'}",
         f"(b) test-like board gain {p['test_like_board']['mean']:+.3f}: {'holds' if p['clauses']['b'] else 'fails'}",
         f"(c) arm minus anti {p['arm_minus_anti']['mean']:+.3f} ({p['arm_minus_anti']['ahead']} of 10 ahead): "
         f"{'holds' if p['clauses']['c'] else 'fails'}",
         f"(d) cell_state gain {p['cell_state']['mean']:+.4f}: {'holds' if p['clauses']['d'] else 'fails'}",
         f"(e) merged gain without e725 {p['without_e725']['mean']:+.3f}: {'holds' if p['clauses']['e'] else 'fails'}",
         f"Outcome: {res['outcome']}. Counted: {', '.join(res['counted'])} (J5 diagnostics).",
         "Reported, deciding nothing (per task: board gain over C0; Kish pool use C0 / arm and library ratio per "
         "bracket):"]
    for task, per in res["reported"].items():
        for n, v in per.items():
            b, a = v["brackets"]["before"], v["brackets"]["after"]
            L.append(f"  {task} {n}: {v['board_gain']:+.3f} ({v['ahead']} ahead); before "
                     f"{_fmt(b['kish_frac_c0'], '.3f')} / {_fmt(b['kish_frac_arm'], '.3f')}, "
                     f"{_fmt(b['library_ratio'], '.3f')}; after "
                     f"{_fmt(a['kish_frac_c0'], '.3f')} / {_fmt(a['kish_frac_arm'], '.3f')}, "
                     f"{_fmt(a['library_ratio'], '.3f')}")
    L.append("  e75 and e725 apart: " + "; ".join(f"{k} arm {v['arm']:+.3f}, anti {v['anti']:+.3f}"
                                                   for k, v in res["e75_e725_apart"].items()))
    L.append("Type shares drawn against C0's per bracket are in reading.json (reported). Nothing here shows E7.75, the "
             "live scorer's build or one specimen's run effects.")
    return "\n".join(L) + "\n"


def load_runs(runs: Path) -> dict[str, list[dict]]:
    out = {}
    for name in nt.all_tasks():
        folder = Path(runs) / "run" / name
        recs = [json.loads(p.read_text()) for p in sorted(folder.glob("seed*.json"))] if folder.exists() else []
        out[name] = sorted(recs, key=lambda r: r["seed"])
    return out


def provenance() -> dict:
    s = repo_state(REPO)
    return {"head": s.commit, "clean": s.unclean == (), "script_sha256": _sha_file(Path(__file__)),
            "t2_mmflow_sha256": _sha_file(Path(mf.__file__)), "harness_sha256": _sha_file(Path(harness.__file__))}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    if not has_row():
        sys.exit("E54 has no row in experiments/harness/TRIES.md")
    runs = Path(a.runs)
    gate, why = load_gate(out, a.commit)
    why = [*why, *repro_problems(runs, a.commit), *commit_problems(a.commit), *exclusion_problems(),
           *fit_problems(runs, a.commit, task_folds().values())]
    res = {**reading(load_runs(runs), gate, why, a.commit), "reading_code": provenance()}
    _write_new(files[0], json.dumps(res, indent=1, default=float))
    _write_new(files[1], text(res))
    print(f"wrote {files[1]}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--fold", required=True)
    f.add_argument("--device", default="cuda")
    for name in ("repro", "run"):
        p = sub.add_parser(name)
        p.add_argument("--task", required=True, help="a short name of nine_tasks.all_tasks()")
    sub.add_parser("gate")
    sub.add_parser("read")
    for p in sub.choices.values():
        p.add_argument("--commit", required=True, help="the run commit (full hash)")
        p.add_argument("--runs", default=str(RUNS))
        p.add_argument("--data", default=str(harness.DATA_DIR))
        p.add_argument("--out", default=str(RESULTS))
    a = ap.parse_args(argv)
    {"fit": cmd_fit, "repro": cmd_repro, "gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
