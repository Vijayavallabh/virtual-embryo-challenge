#!/usr/bin/env python
"""Try E12 (experiments/harness/TRIES.md): the pre-registered readout of the embryo s2 candidate fgw_shape0_aligned.

The candidate is fgw_shape0 (s2, E7.5's geometry at t_pos = 0) with the moves taken from try E10's aligned coupling, so
its mapped later bracket lands on the side position says (FGW's partners are mirrored in every test draw, row N1). This
script applies the row's rule, fixed before the run, to one harness run of floor, fgw, fgw_shape0 and fgw_shape0_aligned
on the two embryo pseudo-tasks at one clean commit (experiments/harness/run.py run; --runs is its output folder) and to
the test build's record (--test-build-json: experiments/builds/build_test.py embryo --method fgw_shape0 --move-coupling
aligned, which exits non-zero with the cloud guard's flags if the guard refuses the cloud, leaving no record).

The kill list (TRIES row E12): stop with no change proposed if any of these holds.
  (a) the expected better of fgw and the candidate, mean over the two tasks, is below 61.02 (fgw_shape0's 60.87 + 0.15);
  (b) on E7.5 the candidate fails the adoption rule's form against fgw_shape0: paired mean gain above +0.3 and ahead in
      at least 8 of the 10 seeds (vec.adoption.MIN_GAIN, MIN_AHEAD);
  (c) de_score, de_direction, mmd_u or variogram differ from fgw_shape0's beyond rounding (the same cells, so expected
      identical): a skill difference above 1e-6 on any seed of either task;
  (d) neighborhood_mmd skill falls by more than 0.01 on E7.5 against fgw_shape0, mean over the seeds (expected sign:
      down; E11 measured -0.010 with aligned moves);
  (e) the test build's cloud guard flags the cloud: no build record, or one that is not this candidate's (entry
      fgw_shape0 with move_coupling aligned at t_pos 0, on the test board, not a dry run). It replaces the designer's
      pre-check kill (partner cross-correlation below 0.6), which the sharpener struck.
It also reports, against the floor and both entries: each task's mean score and SD, the paired gains and seeds ahead,
the expected better-of-two, the per-metric skill changes, and E7.25's mirrored seeds 3, 4, 5, 8, 9 (partners mirrored in
those draws, row E11) separately from the other five, as the mechanism read, not judged. Neither task counts as
adoption-rule evidence for the candidate (Candidate.tuned_on; both were selected on when fgw_shape0 was chosen), so the
verdict decides nothing by itself: an adoption goes through the captain's own ruling only (own commit, new sha256 in
both places of docs/S2_RULES.md, "Changes after the freeze", by 24 Nov).

Preconditions, checked before any reading: one clean commit across both tasks (vec.harness.check_results), the four
candidates on seeds 0-9, and fgw and fgw_shape0 reproducing the committed emb_shape_v1 run per seed (within 1e-6): the
bar 61.02 was read from that run, so a run that does not reproduce it is held, not read. Every threshold is the row's.

It reads aggregate JSONs only (scores and skills per seed, the build record's method fields), prints and writes no cell
count, and refuses to overwrite earlier outputs (a rerun is a new TRIES row). Run from the repo root, BLAS unset (no
linear algebra here):
  uv run --frozen --no-sync python -I -B experiments/harness/e12_readout.py --runs runs/harness/emb_shape_aligned_v1 \
    --test-build-json runs/builds/e12/T2embryo_s2_aligned_e12.build.json --embryo-tries-before 9
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from vec import adoption, harness, provenance
from vec.boards import REPO_ROOT

BOARD = "T2:embryo:test_interp"
TASK_E75, TASK_E725 = "T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25"
FILES = {TASK_E75: "T2e_E7.25_E8.0_E7.5.json", TASK_E725: "T2e_E6.75_E7.5_E7.25.json"}  # as run.py names them
FLOOR, S1, S2, CAND = harness.FLOOR, "fgw", "fgw_shape0", "fgw_shape0_aligned"
NAMES = (FLOOR, S1, S2, CAND)
REFERENCE = REPO_ROOT / "experiments" / "harness" / "results" / "emb_shape_v1" / "tasks"
OUT = REPO_ROOT / "experiments" / "harness" / "results" / "e12_readout"
OUTPUTS = ("e12_readout.json", "e12_readout.txt")
# fixed in row E12
BAR_BETTER_OF_TWO = 61.02  # fgw_shape0's better-of-two 60.87 + 0.15, mean of the two tasks
ROUNDING = 1e-6  # (c) skill units; same cells give exactly equal metrics
BAR_NEIGHBOURHOOD = -0.01  # (d) neighborhood_mmd skill, mean over seeds on E7.5, against fgw_shape0
REPRODUCE_TOL = 1e-6  # score units, per seed, against emb_shape_v1
SAME_CELL_METRICS = ("de_score", "de_direction", "mmd_u", "variogram")
SHAPE_METRICS = ("d2_shape", "occupancy_dice", "scale_log_ratio", "neighborhood_mmd")
MIRRORED_SEEDS = (3, 4, 5, 8, 9)  # E7.25's seeds whose FGW partners are mirrored (row E11)


def _rows(r: dict) -> list[dict]:
    return sorted(r["rows"], key=lambda row: row["seed"])


def scores(r: dict, name: str) -> np.ndarray:
    return np.array([row["scores"][name] for row in _rows(r)], float)


def skills(r: dict, name: str, metric: str) -> np.ndarray:
    return np.array([row["skills"][name][metric] for row in _rows(r)], float)


def adoption_form(gains: np.ndarray) -> bool:
    """The adoption rule's test for one pseudo-task (vec.adoption.keep_change's 'passes')."""
    return bool(gains.mean() > adoption.MIN_GAIN and int((gains > 0).sum()) >= adoption.MIN_AHEAD)


def better_of_two(runs: dict[str, dict], other: str) -> tuple[dict[str, float], float]:
    """The frozen s2 rule's criterion: for each task the mean over seeds of the larger of fgw's and `other`'s score, and
    the mean over the two tasks."""
    per_task = {k: float(np.maximum(scores(r, S1), scores(r, other)).mean()) for k, r in runs.items()}
    return per_task, float(np.mean(list(per_task.values())))


def build_guard(record: dict | None) -> dict:
    """Kill (e) from the test build's record: the build refuses a flagged cloud and leaves no record, so a missing or
    foreign record is the kill; a record of this candidate says the guard passed."""
    if record is None:
        return {"killed": True, "why": "no test build record: the build refused a flagged cloud or did not complete "
                                       "(its log says which; a crash is not a kill: rerun the build before reading)"}
    rec = record.get("record") or {}
    wanted = {"entry": S2, "board": BOARD, "dry_run": False}
    off = [k for k, v in wanted.items() if record.get(k) != v]
    if off:
        return {"killed": True, "why": f"the test build record is not this candidate's on the test board: {off}"}
    if rec.get("move_coupling") != "aligned" or rec.get("t_pos") != 0.0:
        return {"killed": True, "why": "the test build record is not move_coupling aligned at t_pos 0"}
    return {"killed": False, "why": "the build wrote its record: its cloud guard passed"}


def reproduction(runs: dict[str, dict], reference: dict[str, dict] | None) -> dict:
    """fgw and fgw_shape0 against the committed emb_shape_v1 run, per seed (the row's bar was read from that run)."""
    if reference is None:
        return {"ok": False, "why": "the committed emb_shape_v1 task files were not found"}
    worst = {}
    for k, r in runs.items():
        for name in (S1, S2):
            worst[f"{k} {name}"] = float(np.abs(scores(r, name) - scores(reference[k], name)).max())
    return {"ok": all(v <= REPRODUCE_TOL for v in worst.values()), "max_abs_diff": worst}


def read_kills(runs: dict[str, dict], record: dict | None) -> dict:
    """The row's kill list on the run's per-seed numbers and the test build's record; every key's "killed" says whether
    that criterion stops the proposal."""
    r75, r725 = runs[TASK_E75], runs[TASK_E725]
    bo_c, bo_c_mean = better_of_two(runs, CAND)
    bo_s, bo_s_mean = better_of_two(runs, S2)
    gains = scores(r75, CAND) - scores(r75, S2)
    same = {}
    for m in SAME_CELL_METRICS:
        same[m] = max(float(np.abs(skills(r, CAND, m) - skills(r, S2, m)).max()) for r in (r75, r725))
    nb = float((skills(r75, CAND, "neighborhood_mmd") - skills(r75, S2, "neighborhood_mmd")).mean())
    kills = {
        "a": {"what": "expected better of fgw and the candidate, mean of both tasks", "candidate": bo_c_mean,
              "candidate_per_task": bo_c, "fgw_shape0": bo_s_mean, "fgw_shape0_per_task": bo_s,
              "bar": BAR_BETTER_OF_TWO, "killed": bo_c_mean < BAR_BETTER_OF_TWO},
        "b": {"what": "on E7.5 against fgw_shape0, the adoption rule's form", "mean_gain": float(gains.mean()),
              "seeds_ahead": int((gains > 0).sum()), "need_mean_above": adoption.MIN_GAIN,
              "need_ahead": adoption.MIN_AHEAD, "killed": not adoption_form(gains)},
        "c": {"what": "expression metrics differ from fgw_shape0's beyond rounding (same cells expected)",
              "max_abs_skill_diff": same, "rounding": ROUNDING,
              "killed": any(not v <= ROUNDING for v in same.values())},
        "d": {"what": "neighborhood_mmd skill on E7.5 against fgw_shape0 (expected sign: down; E11 -0.010)",
              "mean_change": nb, "bar": BAR_NEIGHBOURHOOD, "killed": not nb >= BAR_NEIGHBOURHOOD},
        "e": {"what": "the test build's cloud guard", **build_guard(record)},
    }
    return kills


def mechanism(r725: dict) -> dict:
    """E7.25's mirrored seeds 3, 4, 5, 8, 9 separately from the other five (row E11): the candidate's paired gain over
    fgw_shape0 and its shape-skill changes; reported, not judged."""
    seeds = np.array([row["seed"] for row in _rows(r725)])
    gain = scores(r725, CAND) - scores(r725, S2)
    delta = {m: skills(r725, CAND, m) - skills(r725, S2, m) for m in SHAPE_METRICS[:2] + SHAPE_METRICS[3:]}
    ahead_of_fgw = scores(r725, CAND) > scores(r725, S1)

    def part(mask: np.ndarray) -> dict:
        return {"seeds": [int(s) for s in seeds[mask]], "mean_gain_over_fgw_shape0": float(gain[mask].mean()),
                "seeds_ahead_of_fgw_shape0": int((gain[mask] > 0).sum()),
                "seeds_ahead_of_fgw": int(ahead_of_fgw[mask].sum()),
                "mean_skill_change": {m: float(d[mask].mean()) for m, d in delta.items()}}

    mirrored = np.isin(seeds, MIRRORED_SEEDS)
    return {"mirrored_seeds": part(mirrored), "other_seeds": part(~mirrored)}


def task_summary(r: dict) -> dict:
    """One task's report: each candidate's mean and SD; the candidate's paired gain over fgw_shape0, fgw and the floor;
    the per-metric skill means of the candidate and the two entries."""
    out = {"scores": {n: {"mean": float(scores(r, n).mean()), "sd": float(scores(r, n).std(ddof=1))} for n in NAMES}}
    for other in (S2, S1, FLOOR):
        g = scores(r, CAND) - scores(r, other)
        out[f"gain_over_{other}"] = {"mean": float(g.mean()), "seeds_ahead": int((g > 0).sum())}
    out["skills"] = {n: {m: float(skills(r, n, m).mean()) for m in (*SAME_CELL_METRICS, *SHAPE_METRICS)}
                     for n in (S1, S2, CAND)}
    return out


def verdict_line(kills: dict, reproduced: dict) -> str:
    if not reproduced["ok"]:
        return ("HOLD: fgw or fgw_shape0 do not reproduce the committed emb_shape_v1 run per seed, so the row's bar "
                "61.02 was read from another baseline; no reading until this is explained")
    killed = [k for k, v in kills.items() if v["killed"]]
    if killed:
        return f"STOP, no change proposed: kill(s) {', '.join(killed)} (see above)"
    return ("CLEARS THE KILL LIST: it goes to the captain; only his own ruling can replace fgw_shape0 in the frozen "
            "embryo rule (own commit, new sha256 in both places of docs/S2_RULES.md, by 24 Nov)")


def report_text(res: dict, stamp: dict) -> str:
    k = res["kills"]
    lines = ["E12 readout (experiments/harness/TRIES.md, try E12): embryo s2 candidate fgw_shape0_aligned against "
             "fgw_shape0 (s2) and fgw (s1)",
             f"Run {res['run_id']} at commit {res['run_commit'][:7]}, clean {res['run_clean']}, seeds 0-9, both "
             "embryo pseudo-tasks. Try count: "
             f"this run is one counted embryo try ({res['tries_before']} before it, {res['tries_before'] + 1} after).",
             f"Verdict: {res['verdict']}", ""]
    for t, s in res["tasks"].items():
        lines.append(t)
        lines.append("  mean (SD): " + "; ".join(f"{n} {v['mean']:.2f} ({v['sd']:.2f})"
                                                  for n, v in s["scores"].items()))
        lines.append("  candidate's paired gain: " + "; ".join(
            f"over {o} {s[f'gain_over_{o}']['mean']:+.2f}, ahead {s[f'gain_over_{o}']['seeds_ahead']}/10"
            for o in (S2, S1, FLOOR)))
        lines.append("  mean skills (fgw / fgw_shape0 / candidate): " + "; ".join(
            f"{m} {s['skills'][S1][m]:.3f} / {s['skills'][S2][m]:.3f} / {s['skills'][CAND][m]:.3f}"
            for m in (*SHAPE_METRICS, *SAME_CELL_METRICS)))
    a, b, c, d, e = k["a"], k["b"], k["c"], k["d"], k["e"]
    lines += ["", "Kill list (the row's):",
              f"  (a) expected better of fgw and it, mean of both tasks: {a['candidate']:.3f} (bar {a['bar']}; "
              f"fgw_shape0 {a['fgw_shape0']:.3f}); per task "
              + ", ".join(f"{t} {v:.3f}" for t, v in a["candidate_per_task"].items())
              + f" -> {'KILLED' if a['killed'] else 'ok'}",
              f"  (b) on E7.5 against fgw_shape0: paired mean {b['mean_gain']:+.3f} (need above "
              f"+{b['need_mean_above']}), ahead in {b['seeds_ahead']} of 10 (need {b['need_ahead']}) -> "
              f"{'KILLED' if b['killed'] else 'ok'}",
              "  (c) expression skill difference from fgw_shape0, largest on any seed of either task: "
              + ", ".join(f"{m} {v:.2e}" for m, v in c["max_abs_skill_diff"].items())
              + f" (rounding {c['rounding']:g}) -> {'KILLED' if c['killed'] else 'ok'}",
              f"  (d) neighborhood_mmd skill change on E7.5 against fgw_shape0: {d['mean_change']:+.4f} (bar "
              f"{d['bar']}; expected sign down, E11 -0.010) -> {'KILLED' if d['killed'] else 'ok'}",
              f"  (e) test build's cloud guard: {e['why']} -> {'KILLED' if e['killed'] else 'ok'}", "",
              "Mechanism read, E7.25's mirrored seeds 3, 4, 5, 8, 9 against the other five (reported, not judged):"]
    for label, m in res["mechanism_e725"].items():
        lines.append(f"  {label}: seeds {m['seeds']}: gain over fgw_shape0 {m['mean_gain_over_fgw_shape0']:+.3f}, "
                     f"ahead {m['seeds_ahead_of_fgw_shape0']}/{len(m['seeds'])}; ahead of fgw "
                     f"{m['seeds_ahead_of_fgw']}/{len(m['seeds'])}; skill change "
                     + ", ".join(f"{a} {v:+.4f}" for a, v in m["mean_skill_change"].items()))
    rep = res["reproduction"]
    lines += ["", ("Reproduction of emb_shape_v1 (fgw, fgw_shape0, per seed, tolerance "
                   f"{REPRODUCE_TOL:g}): {'ok' if rep['ok'] else 'FAILED'} "
                   + (f"(largest {max(rep['max_abs_diff'].values()):.2e})" if rep.get("max_abs_diff") else
                      f"({rep.get('why')})")),
              f"Commit {stamp['commit']}, tree clean {stamp['tree_clean']}, UTC {stamp['utc']}"]
    return "\n".join(lines) + "\n"


def _load(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", required=True, help="the harness run's output folder (run.py run --out)")
    ap.add_argument("--test-build-json", required=True, help="the test build's .build.json record")
    ap.add_argument("--embryo-tries-before", type=int, required=True,
                    help="the embryo board's counted tries before this one, from TRIES.md's running totals, for the "
                         "report's try count")
    ap.add_argument("--reference", default=str(REFERENCE), help="the committed emb_shape_v1 task files")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} exists and is not empty: E12's readout is not revisited; a rerun needs a new TRIES row")
    paths = [out / n for n in OUTPUTS]
    exempt = [*paths, *(provenance.stamp_path(p) for p in paths)]
    started = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    runs_dir = Path(a.runs)
    absent = [FILES[k] for k in (TASK_E75, TASK_E725) if not (runs_dir / FILES[k]).exists()]
    if absent:
        sys.exit(f"the run lacks {absent}: run.py writes <task>.failed.txt for a pseudo-task that failed (the cloud "
                 "guard refuses a flagged build there too); E12 needs both tasks")
    results = [_load(runs_dir / FILES[k]) for k in (TASK_E75, TASK_E725)]
    if {r["task"]["key"] for r in results} != {TASK_E75, TASK_E725}:
        sys.exit("--runs does not hold the two embryo pseudo-tasks")
    missing = sorted({n for r in results for n in NAMES if n not in r["candidates"]})
    if missing:
        sys.exit(f"the run lacks {missing}: E12 is one run of {NAMES} at one clean commit")
    configured = harness.load_tasks()[BOARD]
    try:
        harness.check_results(results, configured, config_sha256=harness._sha256(harness.CONFIG))
    except harness.InconsistentResults as e:
        sys.exit(f"E12 reads one run of the four candidates at one clean commit: {e}")
    runs = {r["task"]["key"]: r for r in results}
    build_path = Path(a.test_build_json)
    record = _load(build_path) if build_path.exists() else None
    ref_paths = {k: Path(a.reference) / FILES[k] for k in FILES}
    reference = {k: _load(p) for k, p in ref_paths.items()} if all(p.exists() for p in ref_paths.values()) else None
    kills = read_kills(runs, record)
    reproduced = reproduction(runs, reference)
    commit = results[0]["provenance"]["commit"]
    res = {"what": "E12 readout: the row's kill list and report for fgw_shape0_aligned", "board": BOARD,
           "candidate": CAND, "entries": {"s2": S2, "s1": S1}, "run_id": runs_dir.name, "run_commit": commit,
           "run_clean": bool(results[0]["provenance"]["tree_clean"]), "seeds": results[0]["seeds"],
           "tries_before": a.embryo_tries_before, "tries_after": a.embryo_tries_before + 1,
           "rule": {"better_of_two_bar": BAR_BETTER_OF_TWO, "adoption_min_gain": adoption.MIN_GAIN,
                    "adoption_min_ahead": adoption.MIN_AHEAD, "rounding": ROUNDING,
                    "neighbourhood_bar": BAR_NEIGHBOURHOOD, "reproduce_tol": REPRODUCE_TOL},
           "kills": kills, "reproduction": reproduced, "verdict": verdict_line(kills, reproduced),
           "tasks": {k: task_summary(r) for k, r in runs.items()}, "mechanism_e725": mechanism(runs[TASK_E725])}
    inputs = {str(runs_dir / FILES[k]): harness._sha256(runs_dir / FILES[k]) for k in FILES}
    inputs[str(build_path)] = harness._sha256(build_path) if build_path.exists() else None
    inputs |= {str(p): harness._sha256(p) for p in ref_paths.values() if p.exists()}
    ended = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    clean = (None if started[1] is None or ended[1] is None
             else not started[1] and not ended[1] and started[0] == ended[0])
    stamp = {"commit": ended[0], "tree_clean": clean, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    text = report_text(res, stamp)
    out.mkdir(parents=True, exist_ok=True)
    (out / "e12_readout.json").write_text(json.dumps({**res, "inputs": inputs, "provenance": stamp}, indent=1) + "\n")
    (out / "e12_readout.txt").write_text(text)
    for p in paths:
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
