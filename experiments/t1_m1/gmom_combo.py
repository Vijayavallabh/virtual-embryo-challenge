#!/usr/bin/env python
"""g_mom on top of M1: diagnostic N22 (target-free) and the labelled read T26 (experiments/harness/TRIES.md).

g_mom (T11) was judged on B0's base and passed its controls T20-T25, but the test build would put it on M1 with R1 and
R2 (row R3). Nothing has measured the two edits together. Both steps below are fixed in their rows before they run.

  n22    Diagnostic N22, not a try: reads no target and no score, and builds no board file. On the test pair (start
         E9.5, reference E10.5, Feng's step 2/5) and the backtest pair (start E8.5, reference E9.5, Feng's step 1/5), on
         all of the reference's fitted-category cells: g_mom's factors applied to the counts as they are, and to the
         counts after M1's Feng direction (R1's held labels and R2's held genes as the build's). In each category the
         realised log change per metacell and gene is regressed through the origin on the factor asked for. The rule
         compares the two slopes: on cm and endocardium, at both pairs, |with M1 / without M1 - 1| must be at most 0.10.
         A second application to the unmoved counts on another stream reports how much of that ratio is rounding noise.
  run    One E10.5 cut of try T26 (run once per cut). `run` refuses to start unless N22 passed. Arms per seed on the
         harness's draws and scoring: b0 and m1_feng (the gates: experiments/t1_m1/results/diag_v1's per-seed scores,
         within 0.02 points), m1_r (M1 with R1 and R2), m1_r_gmom (m1_r with T11's factors as the last count edit).
  read   T26's reading from both cuts' files.

NOT an adoption result: the harness refuses Feng's E9.5 cells on these pseudo-tasks, both cuts are ones T11 was selected
on, and one fit underlies everything. The read can only say whether the edits interfere (the row's raw-unit bar).

Usage (DGX, repo root, BLAS threads capped):
  uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_combo.py n22 --out n22_additivity
  uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_combo.py run --out t26_m1_gmom \
      --task "T1/E9.5->E10.5:chamber"
  uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_combo.py run --out t26_m1_gmom --task "T1/E9.5->E10.5"
  uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_combo.py read --out t26_m1_gmom
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import sparse

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT, load_boards
from vec.core import mean_cp10k
from vec.evaluate import Stage
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
DIAG_V1 = RESULTS / "diag_v1"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_scope"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_hygiene"))
import diagnose  # noqa: E402
import feng_dirs  # noqa: E402
import g_mom_verify as gv  # noqa: E402
import hygiene  # noqa: E402

TASK_CHAMBER, TASK_WHOLE = "T1/E9.5->E10.5:chamber", "T1/E9.5->E10.5"
CUT_FILES = {TASK_CHAMBER: "T1_E9.5_E10.5_chamber.json", TASK_WHOLE: "T1_E9.5_E10.5.json"}
# T11's raw gains over b0 (mean over ten seeds, decreases; the third-wave panel's fix, in T26's row)
T11_RAW_GAIN = {TASK_CHAMBER: {"mmd_u": 0.0084, "variogram": 0.00045},
                TASK_WHOLE: {"mmd_u": 0.0091, "variogram": 0.00051}}
INTERFERENCE_SHARE = 0.5  # the factors' raw gain on M1 under this share of T11's, on a metric and cut, interferes
REPRO_TOL = 0.02  # points per seed, against diag_v1's committed b0 and m1_feng
ARMS = ("b0", "m1_feng", "m1_r", "m1_r_gmom")
N22_TOL = 0.10  # |slope with M1 / slope without - 1|
N22_CATEGORIES = ("cm", "endocardium")  # the rule's; R3's re-registered scope. Others are reported
N22_MIN_CP10K = 1.0  # a metacell's gene enters the regression if at or above this mean CP10k in both before-states
N22_SEED = 20261110
N22_PAIRS = {"test": ("E10.5", "E9.5", 10.5, 12.5), "backtest": ("E9.5", "E8.5", 9.5, 10.5)}  # ref, start, step days
LABEL_N22 = ("diagnostic N22: do g_mom's factors act the same on counts M1 has moved; target-free: no target, no "
             "score, no board, no built file, no cell count")
LABEL_T26 = ("try T26, a labelled diagnostic, not an adoption result: M1 with R1 and R2 with g_mom's factors on top, "
             "in diag_v1's geometry (Feng's E9.5 cells lie inside these pseudo-tasks' external windows; both cuts are "
             "ones T11 was selected on)")


def held_feng(union: Mapping, preregistered: Sequence[str]) -> list[str]:
    """R2's held set for Feng's direction: the pre-registered list joined with R2's union (`build_test._t1_held`)."""
    return sorted(set(preregistered) | set(union["genes"]))


def moved_by_feng(pop, directions: Mapping, fraction: float, groups: Mapping, held: Sequence[str], seed):
    """The population with M1's direction on the cells of Feng's groups and every other cell as it is (R1), through
    `t1_scope.external_shift` on the stream a build gives it."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    covered = np.array([groups.get(lab) in directions for lab in labels], dtype=bool)
    parts = []
    if covered.any():
        idx = np.flatnonzero(covered)
        moved, _ = t1_scope.external_shift(pop[idx].copy(), directions, fraction, groups, held,
                                           label_key=t1_scope.LABEL_KEY, seed=seed)
        parts.append((idx, moved))
    if not covered.all():
        idx = np.flatnonzero(~covered)
        parts.append((idx, pop[idx].copy()))
    return t1_scope._stack_rows(parts, pop.n_obs)


def _metacell_means(C, idx: np.ndarray, member: np.ndarray, k: int) -> np.ndarray:
    C = sparse.csr_matrix(C)
    n_genes = C.shape[1]
    out = np.zeros((k, n_genes), dtype=np.float64)
    for j in range(k):
        cells = idx[member == j]
        if len(cells):
            out[j] = np.asarray(mean_cp10k(C[cells]), dtype=np.float64).ravel()
    return out


def _slope(x: np.ndarray, y: np.ndarray) -> float:
    den = float((x * x).sum())
    return float((x * y).sum() / den) if den > 0 else float("nan")


def additivity_category(pop, moved, no_edit: np.ndarray, factors: ts.StateFactors, held_factors: set[str],
                        category: str, seed: int = N22_SEED) -> dict:
    """One category: the realised log change of g_mom's factors regressed through the origin on the factor asked for, on
    the counts as they are and on the counts M1 moved, over the same (metacell, gene) entries; and the same on another
    stream for the rounding noise."""
    row, logf, blocks = ts._factor_rows(pop, factors, t1_scope.CATEGORY_KEY, "source_obs")
    offset, k, idx = blocks[category]
    idx = idx[~no_edit[idx]]
    member = row[idx] - offset
    held = np.array([g in held_factors for g in factors.genes], dtype=bool)
    x = logf[offset:offset + k].copy()
    x[:, held] = 0.0
    kw = {"rescale": False, "no_edit": no_edit}
    held_list = sorted(held_factors)
    after_pop = ts.apply_state_factors(pop, factors, held_list, [seed, 5], **kw)[0]
    after_noise = ts.apply_state_factors(pop, factors, held_list, [seed, 6], **kw)[0]
    after_moved = ts.apply_state_factors(moved, factors, held_list, [seed, 5], **kw)[0]
    means = {name: _metacell_means(a.layers["counts"], idx, member, k)
             for name, a in (("pop", pop), ("moved", moved), ("after_pop", after_pop), ("after_noise", after_noise),
                             ("after_moved", after_moved))}
    keep = (x != 0) & (means["pop"] >= N22_MIN_CP10K) & (means["moved"] >= N22_MIN_CP10K)
    pc = ts.PSEUDOCOUNT

    def change(after: str, before: str) -> np.ndarray:
        return np.log((means[after][keep] + pc) / (means[before][keep] + pc))

    without = _slope(x[keep], change("after_pop", "pop"))
    with_m1 = _slope(x[keep], change("after_moved", "moved"))
    noise = _slope(x[keep], change("after_noise", "pop"))
    ratio = with_m1 / without if without and np.isfinite(without) else float("nan")
    return {"slope_without_m1": without, "slope_with_m1": with_m1, "ratio": ratio,
            "noise_ratio": noise / without if without and np.isfinite(without) else float("nan"),
            "entries_in_the_regression": int(keep.sum()),
            "passes": bool(np.isfinite(ratio) and abs(ratio - 1) <= N22_TOL)}


def n22_reading(pairs: Mapping[str, Mapping[str, Mapping]]) -> dict:
    """The rule: every category of N22_CATEGORIES present and passing at every pair."""
    missing = [f"{p}:{c}" for p in N22_PAIRS for c in N22_CATEGORIES if c not in pairs.get(p, {})]
    failing = [f"{p}:{c}" for p in N22_PAIRS for c in N22_CATEGORIES
               if c in pairs.get(p, {}) and not pairs[p][c]["passes"]]
    return {"passes": not missing and not failing, "missing": missing, "failing": failing}


def n22_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"Rule: |slope with M1 / slope without - 1| <= {N22_TOL} on {', '.join(N22_CATEGORIES)} "
             f"at both pairs: {'PASSES' if r['passes'] else 'FAILS'}"
             + (f" (failing {r['failing']}, missing {r['missing']})" if not r["passes"] else "")]
    for pair, rows in res["pairs"].items():
        lines.append(f"{pair} pair: Feng's step {res['fractions'][pair]:g}")
        for c, v in rows.items():
            lines.append(f"    {c}: slope without M1 {v['slope_without_m1']:.4f}, with M1 {v['slope_with_m1']:.4f}, "
                         f"ratio {v['ratio']:.4f} (rounding noise alone {v['noise_ratio']:.4f}), "
                         f"{'in the rule' if c in N22_CATEGORIES else 'reported'}")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def run_n22(data: Path, started: Mapping) -> dict:
    spec = t1_scope.load_spec()
    files = {"E8.5": "T1/E8.5_RNA.h5ad", "E9.5": "T1/E9.5_RNA.h5ad", "E10.5": "T1/E10.5_RNA.h5ad"}
    stage = {k: load_stage(data / f, lineage=False) for k, f in files.items()}
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held = held_feng(union, preregistered)
    pairs, fractions = {}, {}
    for pair, (ref_name, start_name, day_a, day_b) in N22_PAIRS.items():
        ref, start = stage[ref_name], stage[start_name]
        src, start_src = f"{ref_name}_RNA", f"{start_name}_RNA"
        cats = [c for c, s in harness.t1_chamber_shares(ref, src).items() if s > 0]
        factors = harness.g_mom_pair_factors(ref, src, start, start_src, categories=cats)[0]
        pop, _labels, no_edit, _ = gv.n21_pool(ref, src, list(factors.categories), spec)
        directions, _ = feng_dirs.feng_directions(list(pop.var_names))
        groups = feng_dirs.label_groups(pop.obs[t1_scope.LABEL_KEY].astype(str).unique())
        fractions[pair] = feng_dirs.bracket_fraction(day_a, day_b)
        moved = moved_by_feng(pop, directions, fractions[pair], groups, held, [N22_SEED, 2])
        pairs[pair] = {c: additivity_category(pop, moved, no_edit, factors, set(preregistered), c)
                       for c in sorted(factors.categories)
                       if (~no_edit & (pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy() == c)).any()}
    res = {"label": LABEL_N22, "pairs": pairs, "fractions": fractions, "reading": n22_reading(pairs),
           "settings": {"tolerance": N22_TOL, "categories": list(N22_CATEGORIES), "min_cp10k": N22_MIN_CP10K,
                        "seed": N22_SEED, "pseudocount": ts.PSEUDOCOUNT}}
    inputs = [data / f for f in files.values()] + [t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ts.R1_CONFIG,
                                                   REPO_ROOT / union["file"]]
    res["stamp"] = gv.stamp(started, inputs)
    return res


def cmd_n22(a) -> None:
    out = RESULTS / a.out
    gv.refuse_existing(out, ("n22.json", "n22.txt"))
    started = gv.tree_state()
    res = run_n22(Path(a.data), started)
    gv.write_outputs(out, "n22", res, n22_text(res))


# ---------------------------------------------------------------------------------------------------------------------
# Try T26


def gmom_hook(ctx: harness.BuildContext, ref, src: str, start, start_src: str, shares: Mapping[str, float],
              held: Sequence[str], covered_labels: set[str]):
    """T11's fold factors as `t1_scope.build_submission`'s last count edit, as `build_test._t1_g_mom` at test: R1's held
    labels unedited, the candidate's own pre-registered held genes, stream [seed, 5], no rescale, and a check that every
    edited label is one Feng's direction covered. ref, start: classified shallow copies, as `harness._g_mom` does."""
    factors = harness._g_mom_factors(ctx, ref, src, start, start_src, shares)

    def hook(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        no_edit = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
        cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
        edited = set(labels[~no_edit & np.isin(cats, list(factors.categories))])
        if not edited <= set(covered_labels):
            extra = sorted(edited - set(covered_labels))
            raise ValueError(f"the hook would edit labels Feng's direction does not cover: {extra}")
        return ts.apply_state_factors(pop, factors, list(held), [ctx.seed, 5], rescale=False, no_edit=no_edit)[0]

    return hook


def m1_r_arms(ctx: harness.BuildContext, directions: Mapping, groups: Mapping, fraction: float, held_m1: Sequence[str],
              held_factors: Sequence[str]) -> dict[str, Stage]:
    """M1 with R1 (no momentum: cells without a Feng group stay as drawn) and R2 (held_m1), with and without g_mom's
    factors as the last edit, from the same draw and streams."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    start, start_src = harness._shallow(ctx.inputs["start"]), ctx.task.role_stems["start"]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    t1_scope.classify(start, source=start_src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    covered = {lab for lab, g in groups.items() if g in directions}
    n = harness.count_for(ctx.board, ref.n_obs)

    def build(post_cells) -> Stage:
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, held_flat=list(held_m1), ref_source=src,
                                           directions=directions, direction_fraction=fraction, direction_groups=groups,
                                           post_cells=post_cells)
        return harness.finish(ctx, ref, dense(sub.X), None)

    return {"m1_r": build(None),
            "m1_r_gmom": build(gmom_hook(ctx, ref, src, start, start_src, shares, held_factors, covered))}


def n22_passed(path: Path) -> bool:
    return path.exists() and bool(json.loads(path.read_text())["reading"]["passes"])


def cmd_run(a) -> None:
    n22 = RESULTS / a.n22 / "n22.json"
    if not n22_passed(n22):
        sys.exit(f"T26 is refused: {gv._rel(n22)} does not hold a passing N22 (its row: T26 needs N22 passed)")
    out = RESULTS / a.out
    name = re.sub(r"[^A-Za-z0-9.]+", "_", a.task).strip("_") + ".json"
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[a.task]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, _ = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                 (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref_stage.n)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    directions, ext = feng_dirs.feng_directions(list(pool.var_names))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = {"b0": b0(ctx), "m1_feng": diagnose.m1_build(ctx, directions, groups, fraction)[0],
                 **m1_r_arms(ctx, directions, groups, fraction, held_m1, preregistered)}
        arms = gv.score_arms(board, target, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res = {"label": LABEL_T26, "task": t.key, "fraction": fraction, "external": ext, "seeds": list(seeds), "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, REPO_ROOT / union["file"]])}
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


def reproduction(rows: Sequence[Mapping], diag_rows: Sequence[Mapping]) -> dict:
    """Per seed, this run's b0 and m1_feng against diag_v1's committed scores."""
    by_seed = {r["seed"]: r for r in diag_rows}
    worst = {"b0": 0.0, "m1_feng": 0.0}
    missing = [r["seed"] for r in rows if r["seed"] not in by_seed]
    for r in rows:
        if r["seed"] in by_seed:
            for arm in worst:
                worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - by_seed[r["seed"]]["scores"][arm]))
    return {"largest_difference": worst, "missing_seeds": missing,
            "holds": not missing and max(worst.values()) <= REPRO_TOL}


def t26_cut(rows: Sequence[Mapping], key: str) -> dict:
    """One cut: the factors' raw gains on M1 against the bar, and the paired points."""
    def arr(arm: str, m: str) -> np.ndarray:
        return np.array([r["arms"][arm]["metrics"][m] for r in rows])

    def score(arm: str) -> np.ndarray:
        return np.array([r["arms"][arm]["score"] for r in rows])

    def group(arm: str, g: str) -> np.ndarray:
        return np.array([r["arms"][arm]["groups"][g] for r in rows])

    raw = {}
    for m, t11 in T11_RAW_GAIN[key].items():
        gain = float((arr("m1_r", m) - arr("m1_r_gmom", m)).mean())
        raw[m] = {"gain_on_m1": gain, "t11_gain": t11, "bar": t11 * INTERFERENCE_SHARE,
                  "share_of_t11": gain / t11, "below_bar": gain < t11 * INTERFERENCE_SHARE}
    d = score("m1_r_gmom") - score("m1_r")
    e = score("m1_r") - score("b0")
    gr = {g: float((group("m1_r_gmom", g) - group("m1_r", g)).mean()) for g in rows[0]["arms"]["m1_r"]["groups"]}
    return {"raw": raw, "interferes": any(v["below_bar"] for v in raw.values()),
            "gmom_on_m1_points": {"mean": float(d.mean()), "ahead": int((d > 0).sum()), "of": len(d)},
            "gmom_on_m1_points_by_group": gr,
            "m1_r_minus_b0_points": {"mean": float(e.mean()), "ahead": int((e > 0).sum()), "of": len(e)},
            "mean_scores": {a: float(score(a).mean()) for a in ARMS}}


def t26_reading(results: Mapping[str, Mapping], diag: Mapping[str, Sequence[Mapping]]) -> dict:
    missing = [k for k in CUT_FILES if k not in results]
    if missing:
        raise ValueError(f"T26 needs both cuts; missing {missing}")
    gates = {k: reproduction(results[k]["rows"], diag[k]) for k in CUT_FILES}
    cuts = {k: t26_cut(results[k]["rows"], k) for k in CUT_FILES}
    if not all(g["holds"] for g in gates.values()):
        outcome = "VOID"
    else:
        outcome = "INTERFERES" if any(c["interferes"] for c in cuts.values()) else "NO INTERFERENCE"
    return {"outcome": outcome, "reproduction": gates, "cuts": cuts}


def t26_text(res: Mapping) -> str:
    lines = [res["label"], f"Rule: the factors' raw gain on M1 under {INTERFERENCE_SHARE:g} of T11's on mmd_u or "
             f"variogram, on either cut, is interference: {res['reading']['outcome']}"]
    for k in (TASK_CHAMBER, TASK_WHOLE):
        g, c = res["reading"]["reproduction"][k], res["reading"]["cuts"][k]
        lines.append(f"{k}{' (primary)' if k == TASK_CHAMBER else ' (scope differs from the test; reported)'}: "
                     f"reproduction of diag_v1 {'yes' if g['holds'] else 'NO'} (largest difference b0 "
                     f"{g['largest_difference']['b0']:.4f}, m1_feng {g['largest_difference']['m1_feng']:.4f})")
        for m, v in c["raw"].items():
            lines.append(f"    {m}: gain on M1 {v['gain_on_m1']:.5f}, T11's {v['t11_gain']:.5f}, share "
                         f"{v['share_of_t11']:.2f}, bar {v['bar']:.5f}, {'BELOW' if v['below_bar'] else 'above'}")
        p, q = c["gmom_on_m1_points"], c["m1_r_minus_b0_points"]
        by = ", ".join(f"{g_} {v:+.2f}" for g_, v in c["gmom_on_m1_points_by_group"].items())
        lines.append(f"    points, reported: factors on M1 {p['mean']:+.2f} (ahead {p['ahead']} of {p['of']}; by group "
                     f"{by}); M1 with R1 and R2 minus b0 {q['mean']:+.2f} (ahead {q['ahead']} of {q['of']})")
        lines.append("    mean scores: " + ", ".join(f"{a} {v:.2f}" for a, v in c["mean_scores"].items()))
    lines.append("Limit: in-sample labelled diagnostic, one fit; never adoption evidence")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = RESULTS / a.out
    gv.refuse_existing(out, ("t26.json", "t26.txt"))
    started = gv.tree_state()
    results = {k: json.loads((out / "tasks" / f).read_text()) for k, f in CUT_FILES.items()}
    diag = {k: json.loads((DIAG_V1 / f).read_text())["rows"] for k, f in CUT_FILES.items()}
    res = {"label": LABEL_T26, "reading": t26_reading(results, diag),
           "stamps": {k: v["stamp"] for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [out / "tasks" / f for f in CUT_FILES.values()])
    gv.write_outputs(out, "t26", res, t26_text(res))
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T26's reproduction gate failed: the reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("n22", help="diagnostic N22, target-free")
    p.add_argument("--out", required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run", help="one cut of try T26")
    p.add_argument("--task", required=True, choices=list(CUT_FILES))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T26's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"n22": cmd_n22, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
