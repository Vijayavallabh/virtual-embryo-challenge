#!/usr/bin/env python
"""Verification of T1 candidate g_mom (try T11, experiments/harness/TRIES.md) before it goes to the captain: its placebo
controls (tries T20 and T21), its composition reading (try T22) and a target-free look at its factors at the test pair
(diagnostic N21). Three steps, every number fixed in their rows before their runs. None of them changes an upload.

  placebo      Tries T20 and T21. Reads a harness run (experiments/harness/run.py run with floor,b0,g_mom,g_mom_flip,
               g_mom_perm on the two counted E10.5 cuts, the same draws and base as g_mom) and applies the row's rule:
               the momentum reading survives only if g_mom beats EACH placebo by a paired mean gain above 0.3 points
               and in at least 8 of 10 seeds on both cuts (the adoption rule, `vec.harness.decide`, with the placebo
               as the base; a loss larger than 0.3 on either cut fails it too). Before reading anything it requires
               the run's b0 to reproduce try T11's committed per-seed scores within 0.02 points and its g_mom within
               0.05 points on both cuts (g_mom carries the fit's PCA and k-means, whose floats can differ between
               runs; b0 does not); otherwise the reading is void and goes to the captain.
  composition  Try T22. The same comparison of g_mom against b0, both built at the truth's own chamber mix, the
               oracle-base construction of experiments/s2_hedge/payoff.py t1-grid (the chamber cut's categories
               without its surrounding cells, renormalised), which removes the composition mismatch of the counted
               cuts (B0 draws E9.5's chamber mix). b0 at that mix (B0's population,
               `vec.t1_state_transport.b0_population`) and g_mom's factors of the fold on top (the fit on the
               reference's own chamber categories, as g_mom's) with R1's held labels unedited, built by the harness's
               builders; scored by the harness's scoring against the chamber cut (the rule, seeds 0-9) and, on the
               same scoring, the whole E10.5 cut (reported, deciding nothing). Rule: the gain must stay at or above
               +0.3 points in at least 8 of 10 seeds on the chamber cut; otherwise the gain is explained by the known
               composition confound (TRIES rule 2's third condition) and g_mom is not proposed. Before anything is
               decided the b0 arm must reproduce the committed t1-grid's oracle-base row
               (experiments/s2_hedge/results/t1_grid/payoff.json) within 0.02 points in the mean over seeds.
  n21          Diagnostic N21, target-free, the test pair. g_mom's factors fitted exactly as the test would fit them
               (start E9.5, reference E10.5, chamber categories; R1's held labels excluded at application), and the
               backtest's fit (start E8.5, reference E9.5) on the categories in both. Reported per category, deciding
               nothing: the size of the factors, the share of cells and genes touched (E10.5's chamber cells, shares
               only), each edited cell's expected depth (unchanged after application, by the keep-depth constant), and
               the concordance of the factor direction with the backtest's fit (the absolute cosine of the leading
               right singular vectors of the two per-metacell factor matrices, and the Spearman correlation of the
               per-gene sizes), with split-half ceilings (one random halving of each stage's cells per category, fixed
               seed; the same fit on each half). It reads the released E8.5 and E9.5 files, the unlocked E10.5 file,
               the held-flat list and the chamber config; no board, no target, no score, no built file; it prints and
               stores no cell count.

Outputs go to experiments/t1_scope/results/<--out>/ (a JSON and a text file per step), stamped with the commit, the
clean-tree flag (vec.provenance's rule, from the start of the run to the write), the UTC time and the sha256 of every
input. A step refuses a folder that already holds its outputs; a rerun needs a new row.

Usage (DGX, repo root, BLAS at most 16 threads, after the row is committed and the harness run has finished for
`placebo`):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
      uv run --frozen --no-sync python -I -B experiments/t1_scope/g_mom_verify.py placebo \
      --runs runs/harness/t20_t21_g_mom_placebo --out t20_t21_placebo
  ... g_mom_verify.py composition --out t22_composition     (try T22; it reads data/raw, not the harness run)
  ... g_mom_verify.py n21 --out n21_g_mom_test_pair         (diagnostic N21)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

from vec import harness, provenance, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import MAX_LOSS, MIN_AHEAD, MIN_GAIN, SEEDS
from vec.boards import GROUPS, RANKED_METRICS, REPO_ROOT, load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
for _d in (HERE, REPO_ROOT / "experiments" / "s2_hedge"):
    sys.path.insert(0, str(_d))
import o1_readout as o1  # noqa: E402
import payoff  # noqa: E402

CUTS = {"T1/E9.5->E10.5": "T1_E9.5_E10.5.json", "T1/E9.5->E10.5:chamber": "T1_E9.5_E10.5_chamber.json"}
TASK_WHOLE, TASK_CHAMBER = tuple(CUTS)
PLACEBOS = {"g_mom_flip": "T20", "g_mom_perm": "T21"}  # candidate -> try
T11_TASKS = REPO_ROOT / "experiments" / "harness" / "results" / "t11_g_mom" / "tasks"
T1_GRID = REPO_ROOT / "experiments" / "s2_hedge" / "results" / "t1_grid" / "payoff.json"
ORACLE_ROW = "bown_cown"  # the grid's oracle-base prediction and truth: the chamber cut's own mix
REPRO_T11_TOL = {"b0": 0.02, "g_mom": 0.05}  # points per seed, T20 and T21's gate against try T11's committed scores
REPRO_SCORE_TOL = 0.02  # points, the mean over seeds, T22's gate against the committed t1-grid
T22_MIN_GAIN = MIN_GAIN  # 0.3 points, the rule's bar, per seed
T22_MIN_SEEDS = MIN_AHEAD  # 8 of the 10 seeds at or above it
SPLIT_SEED = 20261110  # N21's halves: one generator per stage pair, default_rng([SPLIT_SEED, pair]) (test 0, back 1)
N21_APPLY_STREAM = (SPLIT_SEED, 5)  # the one realised application N21 reports (depth only)
DEPTH_TOL = 1e-9  # the largest |expected depth after the edit / the cell's depth - 1| of an edited cell
LABEL_T22 = ("try T22: g_mom against b0 at the truth's own chamber mix (an oracle base, the construction of "
             "experiments/s2_hedge/payoff.py t1-grid); not an adoption result, it can only veto the proposal")
LABEL_N21 = ("diagnostic N21: g_mom's factors at the test pair (start E9.5, reference E10.5), against the backtest's "
             "(start E8.5, reference E9.5); target-free: no target, no score, no board, no built file, no cell count")


def _rel(path: Path) -> str:
    p = Path(path).resolve()
    return p.relative_to(REPO_ROOT).as_posix() if p.is_relative_to(REPO_ROOT) else str(path)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_state() -> dict:
    """HEAD and whether the tree is clean by vec.provenance's rule (None when git cannot tell), as N12 reads it."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    return {"commit": commit, "clean": None if unclean is None else not unclean}


def stamp(started: Mapping, inputs: Sequence[Path]) -> dict:
    """The stamp of a step's outputs: the commit, whether the tree was clean from the start of the run to the write and
    at one commit (vec.provenance's rule), the UTC time and the sha256 of every input file."""
    ended = tree_state()
    clean = (None if started["clean"] is None or ended["clean"] is None
             else started["clean"] and ended["clean"] and started["commit"] == ended["commit"])
    return {"commit": ended["commit"], "tree_clean": clean, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "inputs": {_rel(p): sha256_file(Path(p)) for p in inputs}}


def refuse_existing(out: Path, names: Sequence[str]) -> None:
    if any((out / n).exists() for n in names):
        sys.exit(f"{out} already holds this step's outputs; a rerun needs a new --out and a new TRIES.md row")


def parse_seeds(text: str) -> tuple[int, ...]:
    out: list[int] = []
    for part in text.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return tuple(out)


def write_outputs(out: Path, stem: str, record: dict, text: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{stem}.json").write_text(json.dumps(record, indent=1, default=float) + "\n")
    (out / f"{stem}.txt").write_text(text)
    print(text, end="")


# ---------------------------------------------------------------------------------------------------------------------
# Tries T20 and T21: the placebo controls


def load_run(runs: Path) -> dict[str, dict]:
    """The harness run's two per-cut results (key -> run_task JSON), with b0, g_mom and both placebos in each."""
    out = {}
    for key, fname in CUTS.items():
        path = runs / fname
        if not path.exists():
            sys.exit(f"{path}: no result for {key}; run experiments/harness/run.py with "
                     "floor,b0,g_mom,g_mom_flip,g_mom_perm")
        result = json.loads(path.read_text())
        missing = [c for c in ("b0", "g_mom", *PLACEBOS) if c not in result["candidates"]]
        if missing or result["task"]["key"] != key:
            sys.exit(f"{path}: not {key} with b0, g_mom and both placebos (missing {missing})")
        out[key] = result
    return out


def reproduce_t11(results: Mapping[str, dict]) -> dict:
    """Per cut, seed and arm (b0, g_mom), the absolute difference of score between this run and try T11's committed
    result; passes if every one is within its tolerance (REPRO_T11_TOL points: the arms are the same code, so a larger
    difference means the run is not the comparison T11 made; g_mom carries the fit's PCA and k-means, whose floats can
    differ between runs, b0 does not)."""
    worst = {arm: 0.0 for arm in REPRO_T11_TOL}
    per_cut = {}
    for key, fname in CUTS.items():
        t11 = json.loads((T11_TASKS / fname).read_text())
        by_seed = {row["seed"]: row["scores"] for row in t11["rows"]}
        diffs = {}
        for row in results[key]["rows"]:
            s = row["seed"]
            if s not in by_seed:
                sys.exit(f"{key}: seed {s} is not in try T11's committed result")
            for arm in REPRO_T11_TOL:
                d = abs(float(row["scores"][arm]) - float(by_seed[s][arm]))
                diffs[arm] = max(diffs.get(arm, 0.0), d)
                worst[arm] = max(worst[arm], d)
        per_cut[key] = diffs
    return {"tolerance_points": dict(REPRO_T11_TOL), "largest_difference": worst, "per_cut": per_cut,
            "reproduced": all(worst[arm] <= REPRO_T11_TOL[arm] for arm in worst)}


def read_placebos(results: Mapping[str, dict], configured, config_sha256: str | None) -> dict:
    """For each placebo the adoption rule with the placebo as base and g_mom as change (`harness.decide`: a paired mean
    gain above 0.3 points, ahead in at least 8 of 10 seeds, on both counted cuts, no cut losing more than 0.3), with the
    margin by scored group per cut; and whether the momentum reading survives: every placebo beaten by that rule."""
    out = {}
    for placebo, label in PLACEBOS.items():
        d = harness.decide(list(results.values()), configured, placebo, "g_mom", config_sha256=config_sha256)
        out[placebo] = {"try": label, "decision": d,
                        "group_margins": {k: o1.group_margins(results[k], base=placebo, change="g_mom") for k in CUTS}}
    return {"placebos": out, "survives": all(p["decision"]["keep"] for p in out.values())}


def placebo_text(record: Mapping) -> str:
    st = record["stamp"]
    lines = [f"Tries T20 and T21: g_mom against its placebos on both E10.5 cuts (rule: paired mean gain above "
             f"{MIN_GAIN} points and ahead in at least {MIN_AHEAD} of 10 seeds on both cuts, against each placebo; no "
             f"cut losing more than {MAX_LOSS}): momentum reading "
             f"{'SURVIVES' if record['survives'] else 'DOES NOT SURVIVE'}"]
    rep = record["reproduction_of_t11"]
    lines.append(f"Reproduction of try T11's committed scores (per seed; b0 within {rep['tolerance_points']['b0']}, "
                 f"g_mom within {rep['tolerance_points']['g_mom']} points): "
                 f"{'yes' if rep['reproduced'] else 'NO, the reading is void'} (largest difference b0 "
                 f"{rep['largest_difference']['b0']:.4f}, g_mom {rep['largest_difference']['g_mom']:.4f})")
    for name, p in record["placebos"].items():
        d = p["decision"]
        lines.append(f"{p['try']} {name}: {'beaten' if d['keep'] else 'NOT beaten'} ({d['reason']})")
        for key, v in d["tasks"].items():
            if not v.get("complete"):
                lines.append(f"    {key}: incomplete {v}")
                continue
            m = p["group_margins"][key]
            lines.append(f"    {key}: g_mom minus {name}: paired mean gain {v['mean_gain']:+.2f}, ahead in "
                         f"{v['seeds_ahead']} of 10, passes {v['passes']}; by group "
                         + ", ".join(f"{g} {x:+.2f}" for g, x in m["margins"].items()))
    lines.append(f"Decided on {record['decided_on']}; run {record['run']}; commit {st['commit']}, tree clean "
                 f"{st['tree_clean']}, UTC {st['utc']}")
    return "\n".join(lines) + "\n"


def cmd_placebo(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("placebo.json", "placebo.txt"))
    started = tree_state()
    runs = Path(a.runs)
    results = load_run(runs)
    repro = reproduce_t11(results)
    configured = harness.load_tasks()["T1:test"]
    read = read_placebos(results, configured, harness._sha256(harness.CONFIG))
    if not repro["reproduced"]:
        read["survives"] = False
    inputs = [runs / f for f in CUTS.values()] + [T11_TASKS / f for f in CUTS.values()]
    record = {"label": "tries T20 and T21: g_mom against its placebos on both E10.5 cuts", "run": str(runs),
              "rule": {"min_mean_gain": MIN_GAIN, "min_seeds_ahead": MIN_AHEAD, "max_cut_loss": MAX_LOSS,
                       "every_placebo_on_both_cuts": True},
              "reproduction_of_t11": repro, **read,
              "decided_on": next(iter(read["placebos"].values()))["decision"]["decided_on"],
              "stamp": stamp(started, inputs)}
    write_outputs(out, "placebo", record, placebo_text(record))
    if not repro["reproduced"]:
        sys.exit("the run does not reproduce try T11's committed b0 and g_mom: the reading is void; it goes to the "
                 "captain")


# ---------------------------------------------------------------------------------------------------------------------
# Try T22: the composition reading (the oracle base of the t1-grid)


def oracle_shares(cut_categories: Sequence[str]) -> dict[str, float]:
    """payoff.py t1-grid's oracle base: the chamber categories of the chamber cut (its categories without its
    surrounding cells), renormalised to 1. Nothing here reads the target's expression, only its categories."""
    mix = t1_scope.composition(list(cut_categories))
    chamber = set(t1_scope.load_spec().in_scope("chamber"))
    kept = {k: v for k, v in mix.items() if k in chamber and v > 0}
    total = sum(kept.values())
    return {k: v / total for k, v in kept.items()}


def oracle_arms(ctx: harness.BuildContext, shares: Mapping[str, float], held: Sequence[str]) -> tuple[Stage, Stage]:
    """b0 and g_mom at the mix `shares`, on one build context of the reference (E9.5 with start E8.5, the task's own
    inputs): B0's population (`ts.b0_population`: payoff.py's b0 at that mix, cell for cell), and on the same population
    g_mom's edit: the fold's finished factors (`harness._g_mom_factors`: fitted once per loaded fold on the reference's
    own chamber categories, E9.5's, as g_mom's) applied through `ts.apply_state_factors` on stream [seed, 5], R1's held
    labels unedited. Returns (b0, g_mom) as the harness would finish them (the board's panel, the checks)."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    start, start_src = harness._shallow(ctx.inputs["start"]), ctx.task.role_stems["start"]
    own = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    pop = ts.b0_population(ref, start, shares, n, ctx.seed, list(held), ref_source=src, start_source=start_src,
                           damping=harness.B0_DAMPING)
    genes = ctx.board.genes
    b0 = harness.finish(ctx, ref, dense(pop[:, genes].X).astype(np.float32), None)
    factors = harness._g_mom_factors(ctx, ref, src, start, start_src, own)
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    no_edit = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    moved, _ = ts.apply_state_factors(pop, factors, list(held), [ctx.seed, 5], rescale=False, no_edit=no_edit)
    return b0, harness.finish(ctx, ref, dense(moved[:, genes].X).astype(np.float32), None)


def floor_anchors(t, board, truth: Stage, ref_stage: Stage, seed: int):
    """The floor's anchors for one truth and scoring seed, with the harness's checks (the floor scores exactly 50)."""
    f = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, seed), truth, ref_stage, seed=seed, anchors="computed",
                    build=harness.BUILD, protocol=harness.PROTOCOL)
    missing = sorted(set(RANKED_METRICS[t.task]) - set(f.anchors))
    if abs(f.score - 50.0) > 1e-9 or f.extras["inverted_anchors"] or missing:
        raise ValueError(f"{t.key} seed {seed}: the floor scores {f.score} (inverted {f.extras['inverted_anchors']}, "
                         f"missing anchors {missing})")
    return f


def score_pair(board, truth: Stage, ref_stage: Stage, seed: int, anchors, b0: Stage, gm: Stage) -> dict:
    """b0 and g_mom against one truth's shared anchors: score and group skills of each."""
    out = {}
    for name, pred in (("b0", b0), ("g_mom", gm)):
        r = score_board(board, pred, truth, ref_stage, seed=seed, anchors=anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        out[name] = {"score": float(r.score), "groups": {g: float(v) for g, v in r.group_skill.items()}}
    return out


def t22_reading(rows: Sequence[Mapping], key: str = "chamber") -> dict:
    """Per-seed gain (g_mom minus b0, board points) on one cut, its mean, the seeds at or above the rule's bar (+0.3)
    and above 0, the rule (at least 8 of 10 at or above the bar, on all ten seeds) and the margin by scored group (mean
    over seeds of 100 x weight x the difference of group skills)."""
    gains = {int(r["seed"]): float(r[key]["g_mom"]["score"] - r[key]["b0"]["score"]) for r in rows}
    groups = [(g, w) for g, w, _ in GROUPS["T1"]]
    margins = {g: float(np.mean([100.0 * w * (r[key]["g_mom"]["groups"][g] - r[key]["b0"]["groups"][g]) for r in rows]))
               for g, w in groups}
    v = np.array(list(gains.values()))
    at_or_above = int((v >= T22_MIN_GAIN).sum())
    return {"gains_per_seed": gains, "mean_gain": float(v.mean()), "seeds_at_or_above_bar": at_or_above,
            "seeds_above_zero": int((v > 0).sum()), "n_seeds": len(v), "group_margins": margins,
            "passes": len(v) == len(SEEDS) and at_or_above >= T22_MIN_SEEDS}


def composition_text(record: Mapping) -> str:
    rule = record["rule"]
    ch, wh = record["chamber"], record["whole"]
    st = record["stamp"]
    rep = record["reproduction_of_oracle_b0"]
    verdict = ("incomplete (not all ten seeds)" if not record["complete_seeds"] else
               ("PASSES" if ch["passes"] else
                "FAILS: the gain is explained by the composition confound (TRIES rule 2's third condition) and g_mom "
                "is not proposed"))
    repro = {True: "yes", False: "NO, the reading is void", None: "not read (not all ten seeds)"}[rep["reproduced"]]
    lines = [f"Try T22: {LABEL_T22}",
             f"Rule (chamber cut): the gain at or above {rule['min_gain']:+.1f} points in at least {rule['min_seeds']} "
             f"of {len(SEEDS)} seeds: {verdict}",
             f"Reproduction of the committed t1-grid's oracle-base b0 ({ORACLE_ROW}) on the chamber cut: fresh mean "
             f"{rep['fresh_mean']:.3f}, committed {rep['committed_mean']:.3f} (tolerance {rep['tolerance_points']}): "
             f"{repro}",
             f"Chamber cut, g_mom minus b0 at the oracle mix: mean {ch['mean_gain']:+.2f}; seeds at or above the bar "
             f"{ch['seeds_at_or_above_bar']} of {ch['n_seeds']}, above 0 {ch['seeds_above_zero']}; per seed "
             + ", ".join(f"{g:+.2f}" for g in ch["gains_per_seed"].values())
             + "; by group " + ", ".join(f"{g} {x:+.2f}" for g, x in ch["group_margins"].items()),
             f"Whole cut (reported, deciding nothing; a chamber-only prediction against the whole E10.5 library keeps "
             f"its non-chamber composition gap): mean {wh['mean_gain']:+.2f}; seeds at or above the bar "
             f"{wh['seeds_at_or_above_bar']} of {wh['n_seeds']}; by group "
             + ", ".join(f"{g} {x:+.2f}" for g, x in wh["group_margins"].items()),
             "Oracle mix (shares of the prediction): "
             + ", ".join(f"{k} {v:.3f}" for k, v in sorted(record["oracle_shares"].items())),
             f"Commit {st['commit']}, tree clean {st['tree_clean']}, UTC {st['utc']}"]
    return "\n".join(lines) + "\n"


def cmd_composition(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("composition.json", "composition.txt"))
    started = tree_state()
    seeds = parse_seeds(a.seeds)
    data = Path(a.data)
    tasks = {t.key: t for t in harness.load_tasks()["T1:test"]}
    t_ch, t_wh = tasks[TASK_CHAMBER], tasks[TASK_WHOLE]
    boards = load_boards()
    genes = boards[t_ch.panel_board].genes
    cut, cut_cat = payoff.chamber_cut(t_ch, data, genes)
    shares = oracle_shares(cut_cat)
    grid = json.loads(T1_GRID.read_text())
    want_mix = grid["pred_mix"][ORACLE_ROW]
    if set(want_mix) != set(shares) or max(abs(want_mix[k] - shares[k]) for k in shares) > 1e-9:
        sys.exit("the oracle mix is not the committed t1-grid's: the construction has changed; the reading is void")
    truth_ch = Stage.from_anndata(cut, genes)
    truth_wh = Stage.from_anndata(read_stage(data / t_wh.target, genes), genes)
    ref_stage = Stage.from_anndata(read_stage(data / t_ch.reference, genes), genes)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t_ch.inputs.items()}
    held = list(t1_scope.held_flat_genes())
    board_ch = harness.task_board(t_ch, boards, truth_ch.n, ref_stage.n)
    board_wh = harness.task_board(t_wh, boards, truth_wh.n, ref_stage.n)
    views = {"chamber": (t_ch, board_ch, truth_ch), "whole": (t_wh, board_wh, truth_wh)}
    rows = []
    for s in seeds:
        t0 = time.time()
        row: dict = {"seed": s}
        for key, (t, board, truth) in views.items():
            view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
            ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                       derived={"t1_held_flat": tuple(held)})
            b0, gm = oracle_arms(ctx, shares, held)
            f = floor_anchors(t, board, truth, ref_stage, s)
            row[key] = {"floor": float(f.score), **score_pair(board, truth, ref_stage, s, f.anchors, b0, gm)}
        row["seconds"] = round(time.time() - t0, 1)
        rows.append(row)
        print(f"T22 seed {s} ({row['seconds']} s): chamber gain "
              f"{row['chamber']['g_mom']['score'] - row['chamber']['b0']['score']:+.2f}, whole gain "
              f"{row['whole']['g_mom']['score'] - row['whole']['b0']['score']:+.2f}", flush=True)
    chamber = t22_reading(rows, "chamber")
    whole = t22_reading(rows, "whole")
    fresh = float(np.mean([r["chamber"]["b0"]["score"] for r in rows]))
    committed = float(grid["mean"][ORACLE_ROW][ORACLE_ROW])
    full = tuple(sorted(seeds)) == tuple(SEEDS)
    reproduced = None if not full else abs(fresh - committed) <= REPRO_SCORE_TOL
    inputs_paths = [data / f for f in sorted({t_ch.target, t_ch.reference, *t_ch.inputs.values()})] + [
        t1_scope.HELD_FLAT_GENES, harness.CONFIG, T1_GRID]
    record = {"label": LABEL_T22, "rule": {"min_gain": T22_MIN_GAIN, "min_seeds": T22_MIN_SEEDS, "cut": TASK_CHAMBER},
              "oracle_shares": shares, "chamber": chamber, "whole": whole, "complete_seeds": full,
              "reproduction_of_oracle_b0": {"fresh_mean": fresh, "committed_mean": committed,
                                            "tolerance_points": REPRO_SCORE_TOL, "reproduced": reproduced},
              "seeds": list(seeds), "rows": rows, "stamp": stamp(started, inputs_paths)}
    write_outputs(out, "composition", record, composition_text(record))
    if reproduced is False:
        sys.exit("the oracle-base b0 does not reproduce the committed t1-grid: the reading is void; it goes to the "
                 "captain")


# ---------------------------------------------------------------------------------------------------------------------
# Diagnostic N21: the factors at the test pair


def free_genes(genes: Sequence[str], held: set[str]) -> np.ndarray:
    return np.array([g not in held for g in genes], dtype=bool)


def factor_size(cf: ts.CategoryFactors, free: np.ndarray) -> dict:
    """The size of a finished factor matrix on the genes not held flat: the mean of |log f| over metacells (cell-
    weighted) and genes, the mean over the entries not zeroed and its 90% quantile, the share of entries not zeroed,
    the share of genes touched in at least one metacell, and the share of the entries at the clip."""
    F = np.abs(cf.logf[:, free])
    w = cf.sizes / cf.sizes.sum()
    nz = F > 0
    return {"mean_abs_log_f": float((w * F.mean(1)).sum()),
            "mean_abs_log_f_nonzero": float(F[nz].mean()) if nz.any() else 0.0,
            "q90_abs_log_f_nonzero": float(np.quantile(F[nz], 0.9)) if nz.any() else 0.0,
            "share_entries_touched": float(nz.mean()), "share_genes_touched": float(nz.any(0).mean()),
            "share_entries_at_clip": float((F >= ts.CLIP_AT - 1e-12).mean()), "n_metacells": int(cf.n_metacells)}


def _columns(sf_genes: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """The positions of the genes `keep` (in that order) in the fit's gene order."""
    index = {g: i for i, g in enumerate(sf_genes)}
    return np.array([index[g] for g in keep], dtype=np.int64)


def leading_axis(cf: ts.CategoryFactors, genes: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """The per-category factor direction: the first right singular vector (gene loadings, unit norm, on the genes `keep`
    in that order) of the metacell factor matrix, each metacell weighted by the square root of its size; zeros when the
    matrix is all zero (no direction)."""
    M = np.sqrt(cf.sizes.astype(np.float64))[:, None] * cf.logf[:, _columns(genes, keep)]
    if not M.any():
        return np.zeros(len(keep))
    return np.linalg.svd(M, full_matrices=False)[2][0]


def axis_cosine(u: np.ndarray, v: np.ndarray) -> float:
    """|cosine| of two gene-space axes on the same genes (the sign of a singular vector is arbitrary); NaN when either
    is zero."""
    a, b = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    return abs(float(np.dot(u, v)) / (a * b)) if a > 0 and b > 0 else float("nan")


def gene_size(cf: ts.CategoryFactors, genes: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """The per-gene size of the factors: the cell-weighted standard deviation over metacells of each kept gene's log
    factor (the genes `keep`, in that order)."""
    F = cf.logf[:, _columns(genes, keep)]
    w = cf.sizes / cf.sizes.sum()
    mean = (w[:, None] * F).sum(0)
    return np.sqrt((w[:, None] * (F - mean) ** 2).sum(0))


def profile_spearman(x: np.ndarray, y: np.ndarray) -> float:
    """The Spearman correlation of two per-gene size profiles (NaN when either is constant)."""
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(spearmanr(x, y)[0])


def concordance(a: ts.StateFactors, b: ts.StateFactors, category: str, held: set[str]) -> dict:
    """The concordance of two fits on one category: on the genes both fits hold and the held-flat list leaves free (the
    same genes in the same order for both), the |cosine| of their leading axes (their absolute value: an axis has no
    sign) and the Spearman correlation of their per-gene sizes over the genes with a non-zero size in either fit (the
    many genes with no factor at all would otherwise tie and drive it)."""
    keep = sorted((set(a.genes) & set(b.genes)) - held)
    ca, cb = a.categories[category], b.categories[category]
    if not keep:
        return {"axis_abs_cosine": float("nan"), "size_spearman": float("nan")}
    sa, sb = gene_size(ca, a.genes, keep), gene_size(cb, b.genes, keep)
    moved = (sa > 0) | (sb > 0)
    spearman = profile_spearman(sa[moved], sb[moved]) if moved.sum() >= 3 else float("nan")
    return {"axis_abs_cosine": axis_cosine(leading_axis(ca, a.genes, keep), leading_axis(cb, b.genes, keep)),
            "size_spearman": spearman}


def halves(adata, rng: np.random.Generator) -> np.ndarray:
    """A random half of each category's cells (boolean mask over the rows of a classified stage; floor(n / 2) cells of
    each category in half A), the categories taken in sorted order from one generator."""
    cats = adata.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    mask = np.zeros(adata.n_obs, dtype=bool)
    for c in sorted(set(cats)):
        idx = np.flatnonzero(cats == c)
        mask[rng.permutation(idx)[:len(idx) // 2]] = True
    return mask


def classified(a, source: str, spec):
    c = harness._shallow(a)
    t1_scope.classify(c, source=source, spec=spec)
    return c


def split_half_fits(ref, ref_src: str, start, start_src: str, categories: Sequence[str], pair: int, spec):
    """The same fit on two random halves of each stage's cells of each category: halves drawn from
    default_rng([SPLIT_SEED, pair]), start stage first, then the reference, categories in sorted order. Returns (fit A,
    fit B), each `harness.g_mom_pair_factors` on the cells of its half, for the given categories."""
    rng = np.random.default_rng([SPLIT_SEED, pair])
    ms, mr = halves(classified(start, start_src, spec), rng), halves(classified(ref, ref_src, spec), rng)
    fit_a = harness.g_mom_pair_factors(ref[mr].copy(), ref_src, start[ms].copy(), start_src, categories=categories)[0]
    fit_b = harness.g_mom_pair_factors(ref[~mr].copy(), ref_src, start[~ms].copy(), start_src,
                                       categories=categories)[0]
    return fit_a, fit_b


def expected_depth_ratio(C, factors: np.ndarray, row: np.ndarray, held: np.ndarray) -> np.ndarray:
    """Per cell of counts C, the expected total count after `vec.core.apply_factors` (keep_depth) with these factors,
    over the cell's own total: c * f * s on the free genes (s the keep-depth constant), c on the held genes."""
    C = sparse.csr_matrix(C, dtype=np.float64)
    n = C.shape[0]
    row_of = np.repeat(np.arange(n), np.diff(C.indptr))
    f = factors[row[row_of], C.indices]
    free = ~held[C.indices]
    num = np.bincount(row_of, weights=C.data * free, minlength=n)
    den = np.bincount(row_of, weights=C.data * f * free, minlength=n)
    scale = np.divide(num, den, out=np.ones(n), where=den > 0)
    expected = np.bincount(row_of, weights=C.data * np.where(free, f * scale[row_of], 1.0), minlength=n)
    depth = np.asarray(C.sum(1), dtype=np.float64).ravel()
    return np.divide(expected, depth, out=np.ones(n), where=depth > 0)


def n21_pool(ref, src: str, categories: Sequence[str], spec):
    """E10.5's chamber cells of the fitted categories as a population `apply_state_factors` can read (each cell its own
    source cell, its category and its R1 mask); no draw and no cell number leaves this function."""
    pop = classified(ref, src, spec)
    cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    pop = pop[np.isin(cats, list(categories))].copy()
    pop.obs["source_obs"] = np.asarray(pop.obs_names).astype(str)
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    held_labels = ts.r1_unbracketed_labels(np.unique(labels))
    return pop, labels, np.isin(labels, sorted(held_labels)), sorted(held_labels)


def n21_category(pop, labels, no_edit, factors: ts.StateFactors, held_genes: set[str],
                 category: str) -> tuple[dict, np.ndarray]:
    """One category of the test pair's fit on E10.5's chamber cells: the factor size, the shares of the category's cells
    touched and held by R1 (shares only) and the labels held; returns that record and the category's edited cells
    (boolean over the population, for the depth check)."""
    cf = factors.categories[category]
    cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    in_cat = cats == category
    src = pop.obs["source_obs"].astype(str).to_numpy()
    has_factor = np.fromiter((s in cf.index for s in src), dtype=bool, count=len(src))
    edited = in_cat & has_factor & ~no_edit
    out = {"factor_size": factor_size(cf, free_genes(factors.genes, held_genes)),
           "share_of_category_cells_touched": float(edited.sum() / max(in_cat.sum(), 1)),
           "share_of_category_cells_held_by_r1": float((in_cat & no_edit).sum() / max(in_cat.sum(), 1)),
           "held_labels_in_category": sorted(set(labels[in_cat & no_edit]))}
    return out, edited


def depth_checks(pop, no_edit, factors: ts.StateFactors, held_genes: set[str], edited_any: np.ndarray) -> dict:
    """Over the edited cells: the expected depth after the edit against the cell's own depth (`expected_depth_ratio`,
    exact by construction; a move larger than DEPTH_TOL stops the run), and the realised depth ratio after one
    application of `ts.apply_state_factors` on N21_APPLY_STREAM (median, 1% and 99% quantiles over the edited cells);
    whether every held cell's depth is identical; the share of the fitted categories' cells edited at all; and the
    expected and realised depth per category (per_category, over each category's edited cells)."""
    C = sparse.csr_matrix(pop.layers["counts"], dtype=np.int64)
    row, logf, _ = ts._factor_rows(pop, factors, t1_scope.CATEGORY_KEY, "source_obs")
    row = row.copy()
    row[no_edit] = logf.shape[0] - 1
    held = np.array([g in held_genes for g in factors.genes], dtype=bool)
    F = np.exp(logf)
    F[:, held] = 1.0
    ratio = expected_depth_ratio(C, F, row, held)
    moved, _ = ts.apply_state_factors(pop, factors, sorted(held_genes), list(N21_APPLY_STREAM), rescale=False,
                                      no_edit=no_edit)
    C_new = sparse.csr_matrix(moved.layers["counts"], dtype=np.int64)
    depth = np.asarray(C.sum(1), dtype=np.float64).ravel()
    new_depth = np.asarray(C_new.sum(1), dtype=np.float64).ravel()
    realised = np.divide(new_depth, depth, out=np.ones_like(depth), where=depth > 0)
    idx = np.flatnonzero(edited_any)
    worst = float(np.abs(ratio[idx] - 1.0).max()) if len(idx) else 0.0
    if worst > DEPTH_TOL:
        sys.exit(f"an edited cell's expected depth moves by {worst:.2e} (allowed {DEPTH_TOL:g}): the edit does not "
                 "keep depth; do not use these factors")

    def spread(rows: np.ndarray) -> dict:
        return {"median": float(np.median(realised[rows])) if len(rows) else None,
                "q01": float(np.quantile(realised[rows], 0.01)) if len(rows) else None,
                "q99": float(np.quantile(realised[rows], 0.99)) if len(rows) else None}

    cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    per_category = {}
    for c in sorted(factors.categories):
        rows = np.flatnonzero(edited_any & (cats == c))
        largest = float(np.abs(ratio[rows] - 1.0).max()) if len(rows) else 0.0
        per_category[c] = {"max_abs_expected_ratio_minus_1": largest, "realised_depth_ratio": spread(rows)}
    return {"max_abs_expected_ratio_minus_1": worst,
            "share_of_fitted_cells_edited": float(edited_any.mean()) if pop.n_obs else 0.0,
            "realised_depth_ratio_edited": spread(idx),
            "realised_depth_identical_when_held": bool(np.array_equal(new_depth[no_edit], depth[no_edit])),
            "per_category": per_category}


def n21_text(res: Mapping) -> str:
    def num(x, f="{:.3f}"):
        return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f.format(x)

    d = res["depth"]
    lines = [f"N21: {LABEL_N21}. Decides nothing; reported to the captain with the ruling.",
             f"R1's held labels (no edit): {', '.join(res['r1_held_labels']) or 'none'}. Cells edited: "
             f"{num(d['share_of_fitted_cells_edited'], '{:.2f}')} of the fitted categories' cells",
             f"Depth, all edited cells: largest |E[depth after] / depth - 1| {d['max_abs_expected_ratio_minus_1']:.2e} "
             f"(kept by construction); realised depth ratio median {num(d['realised_depth_ratio_edited']['median'])}, "
             f"1% to 99% {num(d['realised_depth_ratio_edited']['q01'])} to "
             f"{num(d['realised_depth_ratio_edited']['q99'])}; held cells' depth identical "
             f"{d['realised_depth_identical_when_held']}"]
    for c, r in res["test_pair"].items():
        fs = r["factor_size"]
        lines.append(f"{c}: {fs['n_metacells']} metacells; mean |log f| {num(fs['mean_abs_log_f'])} (not zeroed "
                     f"{num(fs['mean_abs_log_f_nonzero'])}, 90% quantile {num(fs['q90_abs_log_f_nonzero'])}); entries "
                     f"touched {num(fs['share_entries_touched'])}, genes touched {num(fs['share_genes_touched'])}, "
                     f"entries at the clip {num(fs['share_entries_at_clip'])}; cells touched "
                     f"{num(r['share_of_category_cells_touched'], '{:.2f}')} of the category (held by R1 "
                     f"{num(r['share_of_category_cells_held_by_r1'], '{:.2f}')}: "
                     f"{', '.join(r['held_labels_in_category']) or 'none'}); expected depth moves at most "
                     f"{d['per_category'][c]['max_abs_expected_ratio_minus_1']:.2e}, realised depth ratio median "
                     f"{num(d['per_category'][c]['realised_depth_ratio']['median'])}")
        con = r.get("concordance_with_backtest")
        if con is None:
            lines.append("    not in the backtest's fit: no concordance")
            continue
        lines.append(f"    concordance with the backtest's fit: |cos| of the leading axes "
                     f"{num(con['axis_abs_cosine'])}, Spearman of the per-gene sizes {num(con['size_spearman'])}")
        ceil = r["split_half_ceiling"]
        lines.append(f"    split-half ceilings (test pair; backtest pair; their geometric mean): axis |cos| "
                     f"{num(ceil['axis_abs_cosine_test'])}; {num(ceil['axis_abs_cosine_backtest'])}; "
                     f"{num(ceil['axis_geometric_mean'])}; size Spearman {num(ceil['size_spearman_test'])}; "
                     f"{num(ceil['size_spearman_backtest'])}")
    st = res["stamp"]
    lines.append(f"Commit {st['commit']}, tree clean {st['tree_clean']}, UTC {st['utc']}")
    return "\n".join(lines) + "\n"


def run_n21(data: Path, started: Mapping) -> dict:
    spec = t1_scope.load_spec()
    files = {"E8.5": "T1/E8.5_RNA.h5ad", "E9.5": "T1/E9.5_RNA.h5ad", "E10.5": "T1/E10.5_RNA.h5ad"}
    stage = {k: load_stage(data / f, lineage=False) for k, f in files.items()}
    held_genes = set(t1_scope.held_flat_genes())
    test_cats = [c for c, s in harness.t1_chamber_shares(stage["E10.5"], "E10.5_RNA").items() if s > 0]
    back_cats = [c for c, s in harness.t1_chamber_shares(stage["E9.5"], "E9.5_RNA").items() if s > 0]
    test_f = harness.g_mom_pair_factors(stage["E10.5"], "E10.5_RNA", stage["E9.5"], "E9.5_RNA", categories=test_cats)[0]
    back_f = harness.g_mom_pair_factors(stage["E9.5"], "E9.5_RNA", stage["E8.5"], "E8.5_RNA", categories=back_cats)[0]
    both = sorted(set(test_f.categories) & set(back_f.categories))
    test_a, test_b = split_half_fits(stage["E10.5"], "E10.5_RNA", stage["E9.5"], "E9.5_RNA", test_cats, 0, spec)
    back_a, back_b = split_half_fits(stage["E9.5"], "E9.5_RNA", stage["E8.5"], "E8.5_RNA", back_cats, 1, spec)
    pop, labels, no_edit, held_labels = n21_pool(stage["E10.5"], "E10.5_RNA", list(test_f.categories), spec)
    test_rows, edited_any = {}, np.zeros(pop.n_obs, dtype=bool)
    for c in sorted(test_f.categories):
        row, edited = n21_category(pop, labels, no_edit, test_f, held_genes, c)
        edited_any |= edited
        if c in both:
            row["concordance_with_backtest"] = concordance(test_f, back_f, c, held_genes)
            ca = concordance(test_a, test_b, c, held_genes) if {c} <= set(test_a.categories) & set(test_b.categories) \
                else None
            cb = concordance(back_a, back_b, c, held_genes) if {c} <= set(back_a.categories) & set(back_b.categories) \
                else None
            at = None if ca is None else ca["axis_abs_cosine"]
            ab = None if cb is None else cb["axis_abs_cosine"]
            both_axes = at is not None and ab is not None and np.isfinite(at) and np.isfinite(ab)
            row["split_half_ceiling"] = {
                "axis_abs_cosine_test": at, "axis_abs_cosine_backtest": ab,
                "axis_geometric_mean": float(np.sqrt(at * ab)) if both_axes else None,
                "size_spearman_test": None if ca is None else ca["size_spearman"],
                "size_spearman_backtest": None if cb is None else cb["size_spearman"]}
        test_rows[c] = row
    res = {"label": LABEL_N21, "pairs": {"test": ["E9.5", "E10.5"], "backtest": ["E8.5", "E9.5"]},
           "categories_test_pair": sorted(test_f.categories), "categories_backtest_pair": sorted(back_f.categories),
           "categories_in_both": both, "skipped_test_pair": sorted(test_f.skipped),
           "skipped_backtest_pair": sorted(back_f.skipped), "r1_held_labels": held_labels,
           "test_pair": test_rows, "depth": depth_checks(pop, no_edit, test_f, held_genes, edited_any),
           "settings": {"damping": harness.B0_DAMPING, "split_seed": SPLIT_SEED,
                        "apply_stream": list(N21_APPLY_STREAM), "zero_below": ts.ZERO_BELOW, "clip_at": ts.CLIP_AT}}
    inputs = [data / f for f in files.values()] + [t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ts.R1_CONFIG]
    res["stamp"] = stamp(started, inputs)
    return res


def cmd_n21(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("n21.json", "n21.txt"))
    started = tree_state()
    res = run_n21(Path(a.data), started)
    write_outputs(out, "n21", res, n21_text(res))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("placebo", help="tries T20 and T21, from a harness run's results")
    p.add_argument("--runs", required=True, help="the harness run folder (experiments/harness/run.py run --out)")
    p.add_argument("--out", required=True, help="folder name under experiments/t1_scope/results/")
    p = sub.add_parser("composition", help="try T22: g_mom against b0 at the truth's own chamber mix")
    p.add_argument("--seeds", default="0-9")
    p.add_argument("--out", required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("n21", help="diagnostic N21: the factors at the test pair (target-free)")
    p.add_argument("--out", required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    {"placebo": cmd_placebo, "composition": cmd_composition, "n21": cmd_n21}[a.cmd](a)


if __name__ == "__main__":
    main()
