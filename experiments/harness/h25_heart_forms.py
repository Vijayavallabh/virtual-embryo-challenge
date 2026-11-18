#!/usr/bin/env python
"""Try H25 of experiments/harness/TRIES.md: the two uploaded forms of the heart file, run end to end on the two
heart backtests, with a short count ladder.

No harness run ever applied the heart's size step or the halving (`vec.harness._b0_t2`: "the heart size cannot be
backtested"), so every practice shape figure of the heart was a copy's (docs/PLAN.md decision log 2026-11-18 (k)).
This script builds each uploaded form as the builder does (experiments/builds/build_test.py, `build_heart`: a
harness build, then `vec.t2_build`'s size step, then the checked builder; then experiments/builds/halve.py's rows) on
the two counted heart backtests, and scores it through the harness's own path (`vec.harness.run_task`: bake
protocol, scorer build main, 10 seeds, build seed BUILD_SEED_OFFSET + seed, the floor checked at 50).

Arms (the first two are references and are not counted):
  copy          the copy baseline, `vec.harness.CANDIDATES["copy"]` (a fresh draw of the reference stage at the count
                rule); its scores must reproduce the committed ones (GUARD_FILES), else the reading is void
  b0_trend      the heart's first-upload expression with no size step, `CANDIDATES["b0_trend"]` (row H11's arm)
  copy_half     (b) the copy, then a uniform random half of its rows
  copy_quarter  (c) the copy, then the half applied twice (a uniform random quarter, nested in copy_half's half)
  iso_half      (d) the ruled second file's form: the copy, the size step in the isotropic raw mode, then the half
  s1form        (e) the first upload's form without the half: b0_trend's expression (`vec.harness._b0_trend`, the
                function the builder calls with --trend, with its own budget rule, B0_DAMPING x the last input gap),
                then the size step in the pre-registered mode xy_grow_z_hold
  s1form_half   (f) s1form, then the half

The size step. The size procedure (experiments/coords/size_procedure.py) gives no figure for E10.5 from the stages
before it: its pre-registered settings use only rate intervals starting at or after E9.0 (heart_min_rate_start_t),
and the backtests' inputs have at most one stage there (E9.5), so every rule is inapplicable and the result is
valid=False (`oracle_note`, checked by the tests on the stage days alone). No size model is fitted here. As
experiments/s2_hedge/payoff.py did, the size is an ORACLE: the target stage's full-cloud RMS radius, measured as the
size procedure measures a stage (`size_procedure.measure_stage`: duplicated records dropped). It enters an in-memory
size result (`oracle_size_result`) that `vec.t2_build.apply_size_procedure` applies, the function `size_or_hold` calls
for a valid heart JSON, with hypothesis H1 and match "full", in the recorded mode (xy_grow_z_hold) or the raw mode.
Both modes reach the same 3-D RMS radius, so the two forms differ in their shape, not their size. Every number that
depends on the size (shape and growth, scale_log_ratio, the board score of a sized arm) is labelled "oracle size".
payoff.py scaled by the core's in-plane growth (its rows differed in scope); here the rows differ in z mode, and both
modes of the size procedure aim at the full cloud's 3-D RMS, so that is the quantity the oracle supplies.

The halving. `halve.half_rows(n, seed)`, the function halve.py applies, with seed HALF_SEED_BASE + seed for the half
and QUARTER_SEED_BASE + seed for the second halving of the quarter. The upload's seed, 1000001 (halve.DEFAULT_SEED),
is not reused: it equals the build seed of scoring seed 1, and numpy's default_rng(1000001) is the stream
`sample_cells` draws that seed's rows with, so the half would be correlated with the draw. Both bases lie apart from
every build seed and every scoring stream, and above halve.py's floor (BUILD_SEED_OFFSET). The quarter is the half
applied twice rather than one draw of n // 4, so that it reuses halve.py's own operation and is nested in the half:
copy -> copy_half -> copy_quarter then removes cells and changes nothing else.

The reading (`reading`), fixed in row H25 before the run. Void if a task, seed or arm is missing, the two results come
from different commits or a dirty tree, the copy does not reproduce its committed scores (ten-seed mean within
GUARD_MEAN board points, every seed within GUARD_SEED), or an arm that shares its cells and expression with another
differs from it on an expression-only skill (iso_half from copy_half, s1form from b0_trend). Then, per arm, backtest
and metric, the arm minus the copy paired by seed, and the board scores; and three statements:
  (1) "The ruled form is at the copy's level on shape": iso_half's d2_shape and occupancy_dice skills each within
      SHAPE_TOL of copy_half's (ten-seed means) on both backtests.
  (2) "Halving costs the neighbourhood term": copy_half below copy on neighborhood_mmd skill by more than NBHD_BAR
      (paired mean) with at least MIN_SEEDS of 10 seeds behind, on both backtests. If true, the AI assistant puts one
      priced question to the captain on the second heart file's cell number, with the ladder (copy, copy_half,
      copy_quarter) as its table; if false, nothing follows.
  (3) "The first upload's form falls below a copy in practice": s1form_half below copy on the board score by more
      than SCORE_BAR (paired mean) with at least MIN_SEEDS of 10 seeds behind, on both backtests. Reported for the
      method report; it decides nothing.
The two backtests share one truth (E10.5): they are one fold of evidence, and "both" is not two confirmations.

  run   (DGX, repo root) score the arms on both backtests; one JSON per backtest under --out
  read  apply the reading to the JSONs in --runs; write reading.json and reading.txt under --out
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

from vec import harness, t2_build
from vec.adoption import SEEDS
from vec.evaluate import Stage

HERE = Path(__file__).resolve().parents[2]  # this checkout (the experiments beside the script)
sys.path.insert(0, str(HERE / "experiments" / "builds"))
sys.path.insert(0, str(HERE / "experiments" / "coords"))
import halve  # noqa: E402
import size_procedure as sp  # noqa: E402

BOARD = "T2:heart:test_extrap"
TASK_KEYS = ("T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5")
REFERENCES = ("copy", "b0_trend")  # not counted
ARMS = ("copy_half", "copy_quarter", "iso_half", "s1form", "s1form_half")  # counted: one try per arm and backtest
ALL = (harness.FLOOR, *REFERENCES, *ARMS)
# arm -> (expression source, size mode or None, number of halvings)
FORMS = {
    "copy_half": ("copy", None, 1),
    "copy_quarter": ("copy", None, 2),
    "iso_half": ("copy", "raw", 1),
    "s1form": ("b0_trend", "xy_grow_z_hold", 0),
    "s1form_half": ("b0_trend", "xy_grow_z_hold", 1),
}
RECORDED_MODE = "xy_grow_z_hold"  # the heart size file's pre-registered mode (the first upload's)
HYPOTHESIS = "H1"  # s1's scope hypothesis, as build_heart names it
HALF_SEED_BASE = 2 * harness.BUILD_SEED_OFFSET  # the half's seed: HALF_SEED_BASE + scoring seed
QUARTER_SEED_BASE = 3 * harness.BUILD_SEED_OFFSET  # the quarter's second halving: QUARTER_SEED_BASE + scoring seed
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")
EXPRESSION_ONLY = ("de_score", "de_direction", "mmd_u")
SIZE_DEPENDENT = ("d2_shape", "occupancy_dice", "scale_log_ratio")  # with the shape group and the score of a sized arm
SAME_CELLS = (("iso_half", "copy_half"), ("s1form", "b0_trend"))  # (arm, twin): same rows and expression
TOL = 1e-9
GUARD_MEAN, GUARD_SEED = 0.1, 0.5  # the copy against its committed scores, board points
GUARD_RUN = "experiments/harness/results/trend_heart_v1/tasks"  # row H11's run, at 8bb23bf
SHAPE_TOL = 0.02  # statement (1), skill
NBHD_BAR = 0.02  # statement (2), skill
SCORE_BAR = 0.3  # statement (3), board points
MIN_SEEDS = 8
SOURCE = "experiments/harness/h25_heart_forms.py (try H25)"


def safe(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", key).strip("_")


# ------------------------------------------------------------------------------------------------ the size step
def oracle_note(input_days, min_start_t: float | None = None) -> list[str]:
    """The size procedure's applicable heart rules for a backtest's input stage days (its own `applicable_rules` at its
    pre-registered settings). Empty: no rule applies, the procedure's result is valid=False and gives no figure."""
    lim = sp.PREREGISTERED["heart_min_rate_start_t"] if min_start_t is None else min_start_t
    return sp.applicable_rules(sorted(float(d) for d in input_days), sp.PREREGISTERED["heart_rules"], lim)


def oracle_size_result(target_rms: float, ref_name: str, ref_t: float, ref_rms: float, ref_sd_z: float,
                       target_t: float) -> dict:
    """A heart size result in the size procedure's layout whose only forecast is the oracle: the target's own full
    RMS radius, for scope H1, in the recorded mode (xy_grow_z_hold; its in-plane figure holds z at the reference
    stage's SD, as the procedure does) and in the raw mode. `vec.t2_build.apply_size_procedure` reads it."""
    r = float(target_rms)
    if not (np.isfinite(r) and r > 0):
        raise ValueError(f"the oracle size must be a positive number, got {target_rms!r}")
    if not r > ref_sd_z:
        raise ValueError("the oracle size is not above the reference's z SD: z cannot be held")
    stage = {"name": ref_name, "t": float(ref_t), "rms_full": float(ref_rms), "sd_z_full": float(ref_sd_z)}
    xy = float(np.sqrt(r ** 2 - float(ref_sd_z) ** 2))
    pre = {"z_mode": RECORDED_MODE, "valid": True, "reference": ref_name, "stages": [stage],
           "scope": {HYPOTHESIS: {"pred_rms_full": r, "pred_rms_xy_full": xy}}}
    raw = {"z_mode": "raw", "valid": True, "reference": ref_name, "stages": [stage],
           "scope": {HYPOTHESIS: {"pred_rms_full": r}}}
    return {"board": "heart", "valid": True, "target_t": float(target_t), "timestamp": "oracle size (try H25)",
            "provenance": {"git_commit": None, "size_procedure_sha256": None, "oracle": True},
            "preregistered": pre, "z_modes": {"raw": raw}}


def size_step(coords: np.ndarray, size: dict, mode: str) -> tuple[np.ndarray, dict]:
    """The builder's size step on a backtest: `apply_size_procedure` on the oracle result, scope H1, match "full";
    the recorded mode when mode is xy_grow_z_hold (as `build_heart` without --z-mode), else the named mode."""
    z_mode = None if mode == RECORDED_MODE else mode
    out, rec = t2_build.apply_size_procedure(coords, size, HYPOTHESIS, z_mode=z_mode, expect_board="heart")
    keep = ("z_mode_used", "z_held", "target_rms_full", "rms_before", "rms_after", "scale", "scale_xy",
            "rms_xy_before", "rms_xy_after", "sd_z_cloud")
    return out, {"label": "oracle size", **{k: rec[k] for k in keep if k in rec}}


# ------------------------------------------------------------------------------------------------ the halving
def kept_rows(n: int, seed: int, halvings: int) -> np.ndarray:
    """The rows kept after `halvings` halvings by halve.half_rows (0: every row) at scoring seed `seed`: the half at
    HALF_SEED_BASE + seed, the second halving at QUARTER_SEED_BASE + seed applied to the half's rows."""
    if halvings not in (0, 1, 2):
        raise ValueError(f"halvings must be 0, 1 or 2, got {halvings!r}")
    rows = np.arange(n)
    for base in (HALF_SEED_BASE, QUARTER_SEED_BASE)[:halvings]:
        s = halve.check_seed(base + int(seed))
        rows = rows[halve.half_rows(len(rows), s)]
    return rows


def form(X: np.ndarray, coords: np.ndarray, size: dict | None, mode: str | None, halvings: int,
         seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    """An arm's file from its harness build's expression and positions: the size step (if any), then the rows kept."""
    rec: dict = {}
    C = np.asarray(coords, dtype=np.float64)
    if mode is not None:
        C, rec = size_step(C, size, mode)
    rows = kept_rows(len(C), seed, halvings)
    return np.asarray(X)[rows], C[rows], rec


# ------------------------------------------------------------------------------------------------ the arms
def make_arm(name: str, size: dict, records: dict) -> harness.Candidate:
    """Arm `name` as a harness candidate for one backtest (its oracle `size`); its size records go to `records`."""
    expr, mode, halvings = FORMS[name]
    base = harness._copy if expr == "copy" else harness._b0_trend

    def build(ctx: harness.BuildContext) -> Stage:
        st = base(ctx)  # the harness build: drawn, moved (b0_trend), guarded and checked
        seed = ctx.seed - harness.BUILD_SEED_OFFSET
        X, C, rec = form(st.X, st.coords, size, mode, halvings, seed)
        sub, rep = t2_build.build_submission(ctx.board, X, C, genes=list(ctx.board.genes))
        if not rep.ok:
            raise ValueError(f"{ctx.task.key} {name}: the built file fails its check: {rep.errors}")
        if rec:
            records.setdefault(name, {})[str(seed)] = rec
        return Stage.from_anndata(sub, ctx.board.genes, need_coords=True)

    ref = harness.CANDIDATES[expr]
    return harness.Candidate(name, frozenset({"T2"}), (), build,
                             {"what": f"try H25: {expr}, size step {mode or 'none'} (oracle size), "
                                      f"{halvings} halving(s) by halve.half_rows", "expression": expr,
                              "size_mode": mode, "halvings": halvings, "size": "oracle" if mode else None},
                             requires=ref.requires if ref.requires is not None else harness._heart_extrap_with_start)


def oracle_for(t: harness.PseudoTask, data_dir: Path) -> tuple[dict, dict]:
    """The oracle size result for backtest t, and its record: the target's and the reference's geometry measured by
    the size procedure's own `measure_stage`, and the procedure's applicable rules on t's input days (none)."""
    g_t, _ = sp.measure_stage(data_dir / t.target, harness.STAGE_DAYS[t.target])
    g_r, _ = sp.measure_stage(data_dir / t.reference, harness.STAGE_DAYS[t.reference])
    size = oracle_size_result(g_t.rms_full, g_r.name, g_r.t, g_r.rms_full, g_r.sd_z_full, g_t.t)
    days = [harness.STAGE_DAYS[f] for f in t.inputs.values()]
    rec = {"label": "oracle size", "target_rms_full": g_t.rms_full, "reference_rms_full": g_r.rms_full,
           "reference_sd_z_full": g_r.sd_z_full,
           "log_target_over_reference": float(np.log(g_t.rms_full / g_r.rms_full)),
           "applicable_rules_at_preregistered_settings": oracle_note(days)}
    return size, rec


def tasks() -> list[harness.PseudoTask]:
    by_key = {t.key: t for t in harness.load_tasks()[BOARD]}
    return [by_key[k] for k in TASK_KEYS]


def cmd_run(a) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    data = Path(a.data)
    for t in tasks():
        path = out / f"{safe(t.key)}.json"
        if path.exists():
            sys.exit(f"{path} exists")
        size, size_rec = oracle_for(t, data)
        if size_rec["applicable_rules_at_preregistered_settings"]:
            sys.exit(f"{t.key}: the size procedure has an applicable rule here; this script uses an oracle only "
                     "where it has none")
        records: dict = {}
        cands = [harness.CANDIDATES[n] for n in REFERENCES] + [make_arm(n, size, records) for n in ARMS]
        missing = [c.name for c in cands if not c.applies(t)]
        if missing:
            sys.exit(f"{t.key}: arms that do not apply: {missing}")
        res = harness.run_task(t, cands, seeds=SEEDS, data_dir=data)
        res.update(task_source=SOURCE, run_id=out.name, oracle_size=size_rec, size_records=records)
        path.write_text(json.dumps(res, indent=1, default=float))
        print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def _rows(r: dict) -> list[dict]:
    return sorted(r["rows"], key=lambda x: x["seed"])


def guard(r: dict, committed: dict) -> list[str]:
    """The copy against its committed scores on one backtest: void reasons, empty when it reproduces."""
    mine = {row["seed"]: row["scores"].get("copy") for row in r["rows"]}
    theirs = {row["seed"]: row["scores"].get("copy") for row in committed["rows"]}
    if any(mine.get(s) is None or theirs.get(s) is None for s in SEEDS):
        return [f"{r['task']['key']}: the copy's scores are missing for a seed here or in the committed run"]
    d = np.array([mine[s] - theirs[s] for s in SEEDS], float)
    why = []
    if abs(d.mean()) > GUARD_MEAN:
        why.append(f"{r['task']['key']}: the copy's ten-seed mean is {d.mean():+.3f} from the committed run's")
    if np.abs(d).max() > GUARD_SEED:
        why.append(f"{r['task']['key']}: a seed's copy score is {np.abs(d).max():.3f} from the committed run's")
    return why


def guard_distances(r: dict, committed: dict) -> dict:
    d = np.array([a["scores"]["copy"] - b["scores"]["copy"] for a, b in zip(_rows(r), _rows(committed), strict=True)])
    return {"mean": float(d.mean()), "largest_seed": float(np.abs(d).max())}


def void_reasons(results: dict[str, dict], committed: dict[str, dict]) -> list[str]:
    """Why the reading may not be read; empty when it may."""
    why = []
    for key in TASK_KEYS:
        r = results.get(key)
        if r is None:
            why.append(f"missing backtest {key}")
            continue
        seeds = sorted(row["seed"] for row in r["rows"])
        if seeds != list(SEEDS):
            why.append(f"{key}: seeds {seeds}, not {list(SEEDS)}")
            continue
        bad = False
        for row in r["rows"]:
            gone = [a for a in ALL if a not in row["scores"]]
            if gone:
                why.append(f"{key} seed {row['seed']}: arms missing {gone}")
                bad = True
                continue
            for arm, twin in SAME_CELLS:
                for m in EXPRESSION_ONLY:
                    if abs(row["skills"][arm][m] - row["skills"][twin][m]) > TOL:
                        why.append(f"{key} seed {row['seed']}: {arm} {m} differs from {twin}'s")
        if not bad:
            if key not in committed:
                why.append(f"{key}: no committed run to check the copy against")
            else:
                why += guard(r, committed[key])
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    return why


def paired(r: dict, arm: str, base: str = "copy") -> dict:
    """arm minus base, paired by seed: the board score, each skill and each group skill (mean, seeds behind, range)."""
    rows = _rows(r)

    def one(vals_a, vals_b) -> dict:
        d = np.asarray(vals_a, float) - np.asarray(vals_b, float)
        return {"mean": float(d.mean()), "behind": int((d < 0).sum()), "ahead": int((d > 0).sum()),
                "lo": float(d.min()), "hi": float(d.max())}

    out = {"score": one([x["scores"][arm] for x in rows], [x["scores"][base] for x in rows])}
    for m in METRICS:
        out[m] = one([x["skills"][arm][m] for x in rows], [x["skills"][base][m] for x in rows])
    for g in GROUPS:
        out[g] = one([x["groups"][arm][g] for x in rows], [x["groups"][base][g] for x in rows])
    return out


def mean_of(r: dict, arm: str, what: str = "score") -> float:
    rows = _rows(r)
    if what == "score":
        return float(np.mean([x["scores"][arm] for x in rows]))
    return float(np.mean([x["skills"][arm][what] for x in rows]))


def statements(results: dict[str, dict]) -> dict:
    """The three statements fixed in row H25, each with the per-backtest numbers that decide it."""
    s1, s2, s3 = {}, {}, {}
    for key in TASK_KEYS:
        r = results[key]
        iso = paired(r, "iso_half", "copy_half")
        s1[key] = {m: iso[m]["mean"] for m in ("d2_shape", "occupancy_dice")}
        half = paired(r, "copy_half")["neighborhood_mmd"]
        s2[key] = {"mean": half["mean"], "behind": half["behind"]}
        sc = paired(r, "s1form_half")["score"]
        s3[key] = {"mean": sc["mean"], "behind": sc["behind"]}
    t1 = all(abs(v) <= SHAPE_TOL for d in s1.values() for v in d.values())
    t2 = all(v["mean"] < -NBHD_BAR and v["behind"] >= MIN_SEEDS for v in s2.values())
    t3 = all(v["mean"] < -SCORE_BAR and v["behind"] >= MIN_SEEDS for v in s3.values())
    return {
        "1": {"statement": "The ruled form is at the copy's level on shape", "true": bool(t1), "per_backtest": s1,
              "rule": f"iso_half minus copy_half within {SHAPE_TOL} on d2_shape and occupancy_dice, both backtests",
              "follows": "nothing (reported)"},
        "2": {"statement": "Halving costs the neighbourhood term", "true": bool(t2), "per_backtest": s2,
              "rule": f"copy_half minus copy on neighborhood_mmd skill below -{NBHD_BAR} with at least {MIN_SEEDS} "
                      "of 10 seeds behind, both backtests",
              "follows": ("one priced question to the captain on the second heart file's cell number, with the "
                          "ladder (copy, copy_half, copy_quarter) as its table, decided from that table alone")
              if t2 else "nothing"},
        "3": {"statement": "The first upload's form falls below a copy in practice", "true": bool(t3),
              "per_backtest": s3,
              "rule": f"s1form_half minus copy on the board score below -{SCORE_BAR} with at least {MIN_SEEDS} of 10 "
                      "seeds behind, both backtests",
              "follows": "reported for the method report; decides nothing; no constant is set from it"},
    }


def ladder(results: dict[str, dict]) -> dict:
    """The count ladder (copy, copy_half, copy_quarter): mean board score and neighbourhood skill per backtest."""
    return {key: {arm: {"score": mean_of(results[key], arm), "neighborhood_mmd": mean_of(results[key], arm,
                                                                                             "neighborhood_mmd")}
                  for arm in ("copy", "copy_half", "copy_quarter")} for key in TASK_KEYS}


def reading(results: dict[str, dict], committed: dict[str, dict]) -> dict:
    why = void_reasons(results, committed)
    if why:
        return {"void": True, "why": why, "statements": None,
                "follows": "nothing: the reading is void (no statement is read)"}
    per = {key: {arm: paired(results[key], arm) for arm in ("b0_trend", *ARMS)} for key in TASK_KEYS}
    means = {key: {arm: {"score": mean_of(results[key], arm),
                         **{m: mean_of(results[key], arm, m) for m in METRICS}} for arm in ALL[1:]}
             for key in TASK_KEYS}
    return {"void": False, "why": [], "commit": next(iter(results.values()))["provenance"]["commit"],
            "guard": {key: guard_distances(results[key], committed[key]) for key in TASK_KEYS},
            "b0_trend_vs_committed": {key: _b0_trend_distance(results[key], committed[key]) for key in TASK_KEYS},
            "oracle_size": {key: results[key].get("oracle_size") for key in TASK_KEYS},
            "statements": statements(results), "ladder": ladder(results), "minus_copy": per, "means": means,
            "labels": {"oracle size": list(SIZE_DEPENDENT) + ["shape_growth", "score of iso_half, s1form, "
                                                                              "s1form_half"]}}


def _b0_trend_distance(r: dict, committed: dict) -> dict | None:
    """Reported only (no void): b0_trend's scores against the committed run's."""
    try:
        d = np.array([a["scores"]["b0_trend"] - b["scores"]["b0_trend"]
                      for a, b in zip(_rows(r), _rows(committed), strict=True)])
    except KeyError:
        return None
    return {"mean": float(d.mean()), "largest_seed": float(np.abs(d).max())}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Try H25 (commit {res['commit']}). Sized arms use an ORACLE SIZE (the target's own RMS radius); every "
             "size-dependent number below is labelled so.",
             "Guard (copy against row H11's committed run): " + "; ".join(
                 f"{k} mean {v['mean']:+.4f}, largest seed {v['largest_seed']:.4f}" for k, v in res["guard"].items())]
    for n, st in res["statements"].items():
        lines.append(f"({n}) {st['statement']}: {'TRUE' if st['true'] else 'FALSE'} ({st['rule']}); follows: "
                     f"{st['follows']}")
        for k, v in st["per_backtest"].items():
            lines.append(f"      {k}: " + ", ".join(f"{a} {b:+.4f}" if isinstance(b, float) else f"{a} {b}"
                                                      for a, b in v.items()))
    for key in TASK_KEYS:
        lines.append(f"{key}: mean board scores " + ", ".join(
            f"{a} {m['score']:.2f}{' (oracle size)' if FORMS.get(a, (None, None))[1] else ''}"
            for a, m in res["means"][key].items()))
        for arm, p in res["minus_copy"][key].items():
            sized = " [oracle size]" if FORMS.get(arm, (None, None))[1] else ""
            lines.append(f"  {arm} - copy{sized}: score {p['score']['mean']:+.2f} (behind {p['score']['behind']} "
                         "of 10); skills " + ", ".join(f"{m} {p[m]['mean']:+.4f}" for m in METRICS))
    return "\n".join(lines) + "\n"


def load_dir(path: Path) -> dict[str, dict]:
    out = {}
    for p in sorted(Path(path).glob("*.json")):
        r = json.loads(p.read_text())
        if isinstance(r, dict) and "rows" in r and r.get("task", {}).get("key") in TASK_KEYS:
            out[r["task"]["key"]] = r
    return out


def cmd_read(a) -> None:
    res = reading(load_dir(Path(a.runs)), load_dir(HERE / GUARD_RUN))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
