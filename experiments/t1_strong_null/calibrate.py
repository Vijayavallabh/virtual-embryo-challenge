#!/usr/bin/env python
"""Row T12 of experiments/harness/TRIES.md (T1; a blind bet taken by the captain's ruling of 2026-11-10 (a),
docs/PLAN.md): the calibration of the strong-null down-tail prior, P(a gene is in the truth's down list | its
reference-level bin, its predicted-change bin), computed with the scorer's own truth lists on the released transitions,
by the rule committed in the row before any run.

Three tables, one per released transition cut:
  table_strong_chamber.json  E9.5 -> E10.5 (the strong-null regime: chance 0.39-0.41, committed de_audit), truth cut to
                             the chamber scope as the harness's task T1/E9.5->E10.5:chamber does. This is the production
                             table: configs/t1_strong_null.yaml pins it (null until the build commit after T12's stop
                             rule) and `build_test.py t1 --strong-null` reads it.
  table_strong_whole.json    the same transition against the whole E10.5 truth: reported, deciding nothing.
  table_weak.json            E8.5 -> E9.5 (the weak-null regime: chance 0.02): the only legal table for T12's backtests
                             of the E9.5 -> E10.5 tasks, which read it as the harness's derived input
                             t1_down_prior_weak.

What is read, and when. The T1 stage files as released: E9.5 and E10.5 (the unlocked validation answer, usable for
preparing P3 predictions under Rules s.10 / A20; docs/DATA_POLICY.md "Final phase, 2026-10-28"), E8.5 and E9.5. E10.5
is used here as a truth (for its down list) only for the production table, which enters only the test build. Feng's
clean statistics (X002, manifest-checked, as build_t1 reads them) and the pinned R2 union enter the strong tables'
predictor. No test board, no E12.5 file, no test target, no score. The script runs once, after T12's row is committed
(it refuses to overwrite its outputs), so nothing is chosen by looking at a table and rerunning.

The rule (every constant in vec.t1_scope and the row). Per transition and scoring seed s of vec.adoption.SEEDS (the
bake protocol's working copies, vec.evaluate.working_copies, exactly as the harness draws them): the truth's down list
is the scorer's `de_genes(truth half, reference copy)` second element (BH 5%, |lfc| >= 0.25). A gene's level bin is its
rank by the FULL reference file's mean log1p CP10k (what a test build knows, the 10% copy being the scorer's draw;
vec.t1_scope.down_prior_bins); its predicted-change bin is its rank by dp = pb(prediction at build seed 1,000,000 + s)
- pb(full reference), the prediction being what supplies "the model's predicted change" at the calibration:
  m1      (strong tables, E9.5 -> E10.5) the test build's own predictor in the calibration's geometry: the reference's
          cells re-weighted to their own chamber mix, Feng's E9.5 -> E14.5 direction on Feng's groups' cells at the
          bracket fraction of the step (0.2), R1 and R2's held genes (`build_test.py t1` as built, with E9.5 in E10.5's
          place). Legal for the production table: Feng's E9.5 and E14.5 lie outside the test's window, and this table is
          never used on a backtest.
  scope   (weak table, E8.5 -> E9.5) the harness's b0_scope predictor: E8.5's cells re-weighted to their own chamber
          mix, no direction (no T1 stage precedes E8.5, and Feng's E9.5 lies inside every backtest's window).
Counts are pooled over the seeds (hits: gene-seeds in the (level, change) cell that were in the down list; genes: all
gene-seeds in the cell) and turned into the table by `DownPriorTable.from_counts` (each cell shrunk towards its level's
marginal by DP_SHRINK gene-seeds). The bins and the shrinkage constant are vec.t1_scope's, fixed in the row.

Outputs, in experiments/t1_strong_null/results/calibration/ (no test-board count, limit or anchor is read or written):
one JSON per table (the table, the pooled counts, per-seed n_up / n_dn, the inputs' sha256 and the stamp),
calibration.txt (a readable summary) and calibration.json (all three with the run's stamp), each with a vec.provenance
build stamp. The
stamp: commit, clean-tree flag (vec.provenance's rule, outputs exempt, from the start of the run to the write), UTC and
the inputs' sha256. It refuses to overwrite existing outputs (a rerun is a new row).

Usage (DGX, repo root; BLAS at most 16 threads; about an hour on 16 cores: per transition cut, one Mann-Whitney
truth-list computation of all genes for each of the ten seeds, plus the predictor builds):
  OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/t1_strong_null/calibrate.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Callable, Sequence
from pathlib import Path

import numpy as np

from vec import harness, provenance, t1_scope
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, working_copies
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "calibration"
for _d in (REPO_ROOT / "experiments" / "t1_m1", REPO_ROOT / "experiments" / "t1_hygiene"):
    sys.path.insert(0, str(_d))
import feng_dirs  # noqa: E402
import hygiene  # noqa: E402

BOARD_KEY = "T1:test"  # the pseudo-tasks' test board in configs/pseudo_tasks.yaml; its limits are never read
TRANSITIONS = ({"name": "strong_chamber", "task": "T1/E9.5->E10.5:chamber", "predictor": "m1",
                "file": "table_strong_chamber.json", "use": "production table (pinned by configs/t1_strong_null.yaml)"},
               {"name": "strong_whole", "task": "T1/E9.5->E10.5", "predictor": "m1",
                "file": "table_strong_whole.json", "use": "reported, deciding nothing"},
               {"name": "weak", "task": "T1/E8.5->E9.5", "predictor": "scope",
                "file": "table_weak.json", "use": "the harness's derived input t1_down_prior_weak (T12's backtests)"})
FRACTION = {9.5: feng_dirs.bracket_fraction(9.5, 10.5)}  # the m1 predictor's bracket fraction at its reference's day


def _empty_table() -> t1_scope.DownPriorTable:
    """A table with the row's bin edges and no values (it only gives `down_prior_bins` its edges)."""
    return t1_scope.DownPriorTable(p=np.zeros((len(t1_scope.DP_LEVEL_EDGES) + 1, len(t1_scope.DP_CHANGE_EDGES) + 2)))


def calibrate_transition(target: Stage, ref: Stage, predict: Callable[[int], np.ndarray],
                         seeds: Sequence[int] = SEEDS) -> dict:
    """The (level, change) counts of one transition. target: the truth stage (the whole target file or its chamber cut)
    and ref: the reference, both as the harness loads them (panel genes, all cells, dense); predict(build_seed): the
    prediction's mean log1p CP10k per panel gene at that build seed. Per seed: the bake protocol's working copies, the
    scorer's down list against the reference copy, the prediction at build seed BUILD_SEED_OFFSET + s, and the bins from
    the FULL reference's pseudobulk and dp = pb(prediction) - pb(full reference). Returns the table and the counts."""
    from common.core_metrics import de_genes  # the scorer's own truth lists (pinned veckit)

    pb_ref = np.asarray(ref.X, dtype=np.float64).mean(0)
    bins = _empty_table()
    hits = np.zeros(bins.p.shape)
    genes = np.zeros(bins.p.shape)
    per_seed = []
    for s in seeds:
        truth, _, ref_w = working_copies(target, ref, int(s), SUBSAMPLE_FRAC, harness.PROTOCOL)
        up, dn, _ = de_genes(truth.X, ref_w.X)
        dp = np.asarray(predict(int(harness.BUILD_SEED_OFFSET + s)), dtype=np.float64) - pb_ref
        level, change = t1_scope.down_prior_bins(pb_ref, dp, bins)
        down = np.zeros(pb_ref.size, dtype=bool)
        down[np.asarray(dn, dtype=np.int64)] = True
        np.add.at(genes, (level, change), 1.0)
        np.add.at(hits, (level[down], change[down]), 1.0)
        per_seed.append({"seed": int(s), "n_up": int(len(up)), "n_dn": int(len(dn))})
    table = t1_scope.DownPriorTable.from_counts(hits, genes)
    return {**table.to_dict(), "shrink": t1_scope.DP_SHRINK, "hits": hits.tolist(), "genes": genes.tolist(),
            "seeds": [int(s) for s in seeds], "per_seed": per_seed,
            "mean_n_up": float(np.mean([r["n_up"] for r in per_seed])),
            "mean_n_dn": float(np.mean([r["n_dn"] for r in per_seed]))}


def _predictor_m1(ref_ad, stem: str, board, fraction: float, held: list[str], directions: dict, groups: dict):
    """The prediction of the production table's predicted-change bin: `build_test.py t1`'s own steps in the
    calibration's geometry (see the module docstring), on the reference's real cells, at the board's cell count."""
    shares = harness.t1_chamber_shares(harness._shallow(ref_ad), stem)
    n = harness.count_for(board, ref_ad.n_obs)

    def predict(seed: int) -> np.ndarray:
        sub, _ = t1_scope.build_submission(ref_ad, board, shares, n, seed, held_flat=held, ref_source=stem,
                                           directions=directions, direction_fraction=fraction,
                                           direction_groups=groups)
        return dense(sub.X).astype(np.float64).mean(0)

    return predict


def _predictor_scope(ref_ad, stem: str, board):
    """The prediction of the weak table's predicted-change bin: the harness's b0_scope step alone."""
    shares = harness.t1_chamber_shares(harness._shallow(ref_ad), stem)
    n = harness.count_for(board, ref_ad.n_obs)

    def predict(seed: int) -> np.ndarray:
        sub, _ = t1_scope.build_submission(ref_ad, board, shares, n, seed, ref_source=stem)
        return dense(sub.X).astype(np.float64).mean(0)

    return predict


def _fmt_table(res: dict) -> list[str]:
    p = np.asarray(res["p"])
    lines = ["  level bin (rank edges " + str(res["level_edges"]) + ", last bin: the rest) x change bin (dp rank edges "
             + str(res["change_edges"]) + ", then the middle, then the top " + str(res["change_top"]) + ")"]
    for i, row in enumerate(p):
        lines.append(f"    L{i}: " + "  ".join(f"{x:.3f}" for x in row))
    return lines


def report_text(summary: dict) -> str:
    s = summary["stamp"]
    lines = [f"T12 calibration of the strong-null down-tail prior; commit {s['commit']}, tree clean {s['tree_clean']}, "
             f"UTC {s['utc']}", ""]
    for name, res in summary["tables"].items():
        lines += [f"{name} ({res['use']}): {res['transition']['task']}, predictor {res['predictor']}, mean n_up "
                  f"{res['mean_n_up']:.1f}, mean n_dn {res['mean_n_dn']:.1f}",
                  *_fmt_table(res), ""]
    lines.append(f"Inputs (sha256): {json.dumps(summary['inputs_sha256'], indent=1)}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("T12's calibration has already written its outputs; a rerun needs a new TRIES row")
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    data = Path(a.data)
    tasks = {t.key: t for t in harness.load_tasks()[BOARD_KEY]}
    boards = load_boards()
    inputs: dict[str, str] = {}
    union = hygiene.load_union()
    held = sorted(set(t1_scope.held_flat_genes()) | set(union["genes"]))
    inputs["union_sha256"] = union["sha256"]
    tables = {}
    for spec in TRANSITIONS:
        t = tasks[spec["task"]]
        genes = boards[t.panel_board].genes
        if t.target_filter is not None:
            target, _ = harness.filter_target(t, data, genes)
        else:
            target = Stage.from_anndata(read_stage(data / t.target, genes), genes)
        ref = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
        board = harness._candidate_board(harness.task_board(t, boards, target.n, ref.n))
        stem = harness.stem(t.reference)
        ref_ad = load_stage(data / t.reference, lineage=False)
        for f in (t.target, t.reference):
            inputs[f] = hygiene.sha256_file(data / f)
        if spec["predictor"] == "m1":
            directions, ext = feng_dirs.feng_directions(list(ref_ad.var_names))
            inputs["feng_stats"] = ext["sha256"]
            groups = feng_dirs.label_groups(ref_ad.obs["celltype"].astype(str).unique())
            predict = _predictor_m1(ref_ad, stem, board, FRACTION[harness.STAGE_DAYS[t.reference]], held, directions,
                                    groups)
            predictor = ("m1: the reference's cells re-weighted to their own chamber mix, Feng's E9.5 -> E14.5 "
                         "direction on Feng's groups' cells at the step's bracket fraction, R1 and R2's held genes")
        else:
            predict = _predictor_scope(ref_ad, stem, board)
            predictor = "scope: the reference's cells re-weighted to their own chamber mix (b0_scope), no shift"
        res = calibrate_transition(target, ref, predict)
        tables[spec["name"]] = {"file": spec["file"], "use": spec["use"], "predictor": predictor,
                                "transition": {"task": t.key, "reference": t.reference, "truth": t.target,
                                               "truth_filter": dict(t.target_filter) if t.target_filter else None},
                                **res}
        print(f"{spec['name']}: {spec['task']} done (mean n_up {res['mean_n_up']:.1f}, "
              f"mean n_dn {res['mean_n_dn']:.1f})", flush=True)
    commit_end, unclean_end = provenance.git_state(provenance.REPO_ROOT)
    stamp = {"commit": commit, "commit_end": commit_end,
             "tree_clean": unclean == [] and unclean_end == [] and commit == commit_end,
             "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    summary = {"row": "T12 (experiments/harness/TRIES.md)", "tables": tables, "inputs_sha256": inputs, "stamp": stamp,
               "constants": {"level_edges": list(t1_scope.DP_LEVEL_EDGES),
                             "change_edges": list(t1_scope.DP_CHANGE_EDGES), "change_top": t1_scope.DP_CHANGE_TOP,
                             "shrink": t1_scope.DP_SHRINK}}
    out.mkdir(parents=True, exist_ok=True)
    for spec in TRANSITIONS:
        path = out / spec["file"]
        path.write_text(json.dumps({**tables[spec["name"]], "inputs_sha256": inputs, "stamp": stamp}, indent=1) + "\n")
        provenance.write_stamp(path)
    (out / "calibration.json").write_text(json.dumps(summary, indent=1) + "\n")
    provenance.write_stamp(out / "calibration.json")
    text = report_text(summary)
    (out / "calibration.txt").write_text(text)
    provenance.write_stamp(out / "calibration.txt")
    print(text, end="")


if __name__ == "__main__":
    main()
