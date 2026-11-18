#!/usr/bin/env python
"""Reading of tries H24 (heart) and E46 (embryo) of experiments/harness/TRIES.md: herded draws on each T2 board's
first-upload entry, scored by `experiments/harness/run.py run` on the board's counted backtests (decision log
2026-11-18 (c) fixes its use before any run).

Why. Try T5 read herded draws on Task 1 (within each category's quota a drawn cell is swapped for an undrawn one when
that brings the sample's mean expression closer to the pool's; `vec.core.herd_rows`). Rows H4 and E9 registered the
same option on the T2 boards (`b0_herd`, `fgw_herd`). On T2 the option is not T5's per-quota form: `sample_cells` and
`interpolate_clouds` herd over every row of the stage's pool (`herd_rows(..., np.arange(N), ...)`), so a swap can
cross types and positions, and herding can move the drawn type shares and positions as well as the mean expression.
This reading puts the option on each board's first-upload entry as the harness builds it: the heart's `b0_trend`
against `b0_trend_herd` (`b0_trend` with the draw herded as `b0_herd`'s), the embryo's `fgw` against `fgw_herd` (the
embryo's first upload is `fgw` itself: `build_test.py embryo --method fgw --expression mixture`, so its herded twin is
the registered `fgw_herd`).

What it reads. Only the harness task files of one board in `--runs` (per-seed scores, skills and metrics), and the
committed task files (`COMMITTED`) that hold the entry's scores, to check that the run reproduces them. It opens no
data file, runs no candidate and no scorer.

Validity, checks that can fail (otherwise the reading is void and decides nothing). (1) The run's task files pass
vec.harness.check_results: one commit, a clean tree, the current config, the bake protocol, seeds 0-9 and their build
seeds; no backtest failed. (2) Every counted backtest of the board (configured, counted, and every arm applies) holds
every arm (the entry, its herded twin and the registered pair: four arms on the heart; two on the embryo, where the
registered pair is the entry and its twin) with a finite score on each of the ten seeds. (3) Where a committed harness
result holds an arm's scores on a backtest (`COMMITTED`; the entry's, and the registered pair's base where it holds
them), the run agrees with them within TOL on every seed; where none holds the entry's, the reading says the
reproduction could not be checked. The harness task files record no machine: the committed runs' machine is named in
the row, and the run must use it (a BLAS or machine difference voids the reading through (3)).

Idle herding. `vec.core.herd_rows` returns the draw unchanged when it repeats rows or nothing is left to swap in (a
board count at or above the pool's rows). A counted backtest on which the herded twin equals its base bit for bit on
every seed (the score and every metric and skill the task file holds) is a backtest on which herding was not applied:
it is reported as idle and its herding is not read there. This test uses no count.

Statement, per board, one of three: vec.adoption.keep_change of the herded entry against the entry on the board's two
counted backtests (both must pass): "the herded entry passes the adoption rule" or "does not pass"; or, where herding
is idle on a counted backtest, "the herded entry's read is not readable" (the rule needs a pass on both backtests, and
an idle backtest reads nothing of herding; the other backtest's numbers are still reported). With each backtest's
paired mean, the paired difference's seed SD, seeds ahead, largest single loss and seeds on which the two arms are
identical. A "does not pass" does not exclude a true gain of about 0.3-0.4 points: the bar (mean > 0.3 and 8 of 10
seeds ahead, on both) has low power there. The registered pair (b0 / b0_herd, fgw / fgw_herd) is reported the same way
and decides nothing. Reported too, descriptive and deciding nothing: the position-only metrics (the shape group:
d2_shape, occupancy_dice, scale_log_ratio) and the shape group's skill of every arm, each with its per-seed paired
difference from the entry (mean and seed SD) beside the entry's own seed-to-seed SD. A herded draw is a different set
of cells, so these differ through the change of draw alone; they show the size of the change, not that herding does or
does not act on position (a scored arm passed the T2 head's cloud guard, which raises otherwise).

  python -I -B experiments/harness/n53_herd_t2_readings.py read --board heart --runs DIR --out DIR
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from vec import harness  # noqa: E402
from vec.adoption import SEEDS, keep_change  # noqa: E402
from vec.boards import GROUPS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

TOL = 1e-9  # per seed, board points: the reproduction of committed scores
SHAPE_GROUP = "shape_growth"
POSITION_ONLY = tuple(m for g, _, ms in GROUPS["T2"] if g == SHAPE_GROUP for m, _ in ms)
OUTCOMES = ("passes", "does not pass", "not readable")


@dataclass(frozen=True)
class Spec:
    test_board: str
    entry: str  # the board's first-upload entry as the harness builds it
    herded: str  # its herded twin
    pair: tuple[str, str]  # the registered pair (rows H4 and E9): base, herded
    committed: Path  # the committed harness task files that hold the entry's scores
    limits: str = ""  # a limit of this board's read

    @property
    def arms(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys((self.entry, self.herded, *self.pair)))

    @property
    def bases(self) -> tuple[str, ...]:
        """The arms whose committed scores the run must reproduce where a committed result holds them (the herded
        arms are new)."""
        return tuple(dict.fromkeys((self.entry, self.pair[0])))


BOARDS = {
    "heart": Spec("T2:heart:test_extrap", "b0_trend", "b0_trend_herd", ("b0", "b0_herd"),
                  HERE / "results" / "trend_heart_v1" / "tasks",  # row H11's run of b0_trend
                  "the heart's backtest T2h/<=E8.75->E10.5 chose B0's composition damping (configs/pseudo_tasks.yaml): "
                  "the entry's own level there is in-sample, the paired difference does not depend on it."),
    "embryo": Spec("T2:embryo:test_interp", "fgw", "fgw_herd", ("fgw", "fgw_herd"),
                   HERE / "results" / "emb_fgw_v1" / "tasks"),  # the embryo entry's run of 2026-11-01
}

LIMITS = ("Limits: herding has no biological content in its intent (it moves the drawn sample's mean panel expression "
          "towards the pool's), but on T2 it swaps cells across the whole pool, not within each type's quota as on "
          "Task 1, so it can move the drawn type shares and positions too; the scorer's own reference working copy "
          "carries much of the de terms' noise, which herding the prediction cannot remove; each board has two counted "
          "backtests and the rule needs a pass on both: a pass is in-sample for nothing but is still two tasks; a "
          "'does not pass' does not exclude a true gain of about 0.3-0.4 points (low power at the bar); the harness "
          "task files record no machine, so the run must use the machine of the committed runs named in the row.")


def load(folder: Path, test_board: str) -> tuple[list[dict], list[str]]:
    """The harness task files of `test_board` in `folder`, and the reasons the run is void found while loading (a
    backtest that failed, a file that is not a task file of the board is skipped)."""
    folder = Path(folder)
    results, why = [], []
    if not folder.is_dir():
        return [], [f"{folder} is not a folder"]
    for p in sorted(folder.glob("*.failed.txt")):
        why.append(f"a backtest failed ({p.name})")
    for p in sorted(folder.glob("*.json")):
        try:
            r = json.loads(p.read_text())
        except (OSError, ValueError):
            why.append(f"{p.name} is not readable JSON")
            continue
        if isinstance(r, dict) and "task" in r and "rows" in r and r["task"].get("test_board") == test_board:
            results.append(r)
    return results, why


def counted(spec: Spec, configured) -> list[str]:
    """The board's counted backtests for this reading: configured, counted (not a diagnostic or a replicate fold),
    where every arm applies and no arm was tuned on it (as vec.harness.decide counts them)."""
    out = []
    for t in configured:
        cands = [harness.CANDIDATES[a] for a in spec.arms]
        if t.counted and all(c.applies(t) for c in cands) and not any(t.key in c.tuned_on for c in cands):
            out.append(t.key)
    return out


def scores_of(r: dict, arm: str) -> dict[int, float] | None:
    """{seed: score} of `arm` in task file `r`, or None if the file does not hold it on every row."""
    if arm not in r.get("candidates", []):
        return None
    out = {}
    for row in r["rows"]:
        if arm not in row.get("scores", {}):
            return None
        out[int(row["seed"])] = float(row["scores"][arm])
    return out


def complete(s: dict[int, float] | None) -> bool:
    return s is not None and sorted(s) == list(SEEDS) and all(math.isfinite(s[k]) for k in SEEDS)


def reproduction(spec: Spec, run: dict[str, dict], keys: list[str], committed_dir: Path) -> tuple[dict, list[str]]:
    """Per counted backtest and base arm (`Spec.bases`), whether the run reproduces the arm's committed scores:
    {"checked": bool, "max_abs": float | None, "file": str | None, "commit": str | None}; and the reasons the run is
    void (a difference above TOL)."""
    committed = {}
    if Path(committed_dir).is_dir():
        for p in sorted(Path(committed_dir).glob("*.json")):
            r = json.loads(p.read_text())
            if isinstance(r, dict) and "task" in r and "rows" in r:
                committed[r["task"]["key"]] = (r, p)
    out, why = {}, []
    for k in keys:
        have = committed.get(k)
        out[k] = {}
        for arm in spec.bases:
            old = scores_of(have[0], arm) if have else None
            if not complete(old):
                out[k][arm] = {"checked": False, "max_abs": None, "file": None, "commit": None}
                continue
            new = scores_of(run[k], arm)
            d = max(abs(new[s] - old[s]) for s in SEEDS)
            out[k][arm] = {"checked": True, "max_abs": d, "file": _rel(have[1]),
                           "commit": have[0].get("provenance", {}).get("commit")}
            if not d <= TOL:
                why.append(f"{k}: {arm} differs from its committed scores ({_rel(have[1])}) by {d:.2e}")
    return out, why


def identical_seeds(r: dict, base: str, change: str) -> int:
    """The number of seeds on which `change` equals `base` bit for bit in task file `r`: the same score and the same
    value of every metric and skill the row holds for either arm."""
    n = 0
    for row in r["rows"]:
        same = row["scores"][base] == row["scores"][change]
        for field in ("metrics", "skills"):
            got = row.get(field) or {}
            if base in got or change in got:
                same = same and got.get(base) == got.get(change)
        n += bool(same)
    return n


def paired(run: dict[str, dict], keys: list[str], base: str, change: str) -> dict:
    """keep_change's verdict of `change` against `base`, with each backtest's seed SD of the paired difference, largest
    single loss (the largest per-seed base - change, 0 if the change is ahead or level on every seed), the seeds on
    which the two arms are identical, and whether herding is idle there (identical on every seed). The outcome is
    "not readable" where it is idle on any counted backtest, else keep_change's."""
    pairs = {k: {s: (scores_of(run[k], base)[s], scores_of(run[k], change)[s]) for s in SEEDS} for k in keys}
    out = keep_change(pairs)
    for k in keys:
        t = out["tasks"][k]
        diffs = np.array([c - b for b, c in (pairs[k][s] for s in SEEDS)])
        t["sd"] = float(diffs.std(ddof=1))
        t["largest_loss"] = max(0.0, float(-diffs.min()))
        t["identical_seeds"] = identical_seeds(run[k], base, change)
        t["idle"] = t["identical_seeds"] == len(SEEDS)
    idle = [k for k in keys if out["tasks"][k]["idle"]]
    outcome = "not readable" if idle else ("passes" if out["keep"] else "does not pass")
    return {"base": base, "change": change, **out, "idle": idle, "outcome": outcome}


def position(run: dict[str, dict], keys: list[str], arms, entry: str) -> dict:
    """Descriptive, per backtest and arm: each position-only metric (raw values) and the shape group's skill (the
    group's within-group weights on the skills), as the mean over seeds, the seed-to-seed SD, and, for an arm other
    than the entry, the mean and seed SD of its per-seed difference from the entry."""
    weights = dict(next(ms for g, _, ms in GROUPS["T2"] if g == SHAPE_GROUP))
    out = {}
    for k in keys:
        rows = sorted(run[k]["rows"], key=lambda row: int(row["seed"]))

        def values(a, rows=rows):
            v = {m: np.array([row["metrics"][a][m] for row in rows], float) for m in POSITION_ONLY}
            v[SHAPE_GROUP] = np.array([sum(w * row["skills"][a][m] for m, w in weights.items()) for row in rows])
            return v

        ve = values(entry)
        out[k] = {}
        for a in arms:
            va = values(a)
            out[k][a] = {m: {"mean": float(va[m].mean()), "sd": float(va[m].std(ddof=1)),
                             **({} if a == entry else {"diff_mean": float((va[m] - ve[m]).mean()),
                                                       "diff_sd": float((va[m] - ve[m]).std(ddof=1))})}
                         for m in (*POSITION_ONLY, SHAPE_GROUP)}
    return out


def reading(board: str, runs: Path, committed_dir: Path | None = None, configured=None) -> dict:
    spec = BOARDS[board]
    committed_dir = spec.committed if committed_dir is None else Path(committed_dir)
    configured = harness.load_tasks()[spec.test_board] if configured is None else configured
    results, why = load(runs, spec.test_board)
    head = {"board": board, "test_board": spec.test_board, "entry": spec.entry, "herded": spec.herded,
            "pair": list(spec.pair), "arms": list(spec.arms), "runs": _rel(runs)}
    if not results:
        return {"void": True, "why": [*why, f"no task file of {spec.test_board} in the run"], **head}
    try:
        harness.check_results(results, configured, harness._sha256(harness.CONFIG))
    except harness.InconsistentResults as e:
        why.append(f"the run's task files do not belong together: {e}")
    keys = counted(spec, configured)
    head["counted"] = keys
    if not keys:
        why.append(f"no counted backtest of {spec.test_board} where every arm applies")
    run = {r["task"]["key"]: r for r in results}
    for k in keys:
        if k not in run:
            why.append(f"{k}: no task file in the run")
            continue
        for a in spec.arms:
            if not complete(scores_of(run[k], a)):
                why.append(f"{k}: {a} is not scored on each of the ten seeds")
    if why:
        return {"void": True, "why": sorted(set(why)), **head}
    rep, why = reproduction(spec, run, keys, committed_dir)
    head["reproduction"] = rep
    if why:
        return {"void": True, "why": sorted(set(why)), **head}
    stmt = paired(run, keys, spec.entry, spec.herded)
    pair = paired(run, keys, *spec.pair)
    commit = results[0]["provenance"]["commit"]
    return {"void": False, "why": [], **head, "commit": commit, "outcome": stmt["outcome"],
            "passes": {"passes": True, "does not pass": False, "not readable": None}[stmt["outcome"]],
            "statement": stmt, "registered_pair": pair,
            "pair_is_the_statement": tuple(spec.pair) == (spec.entry, spec.herded),
            "position_only": position(run, keys, spec.arms, spec.entry)}


def _rows(res: dict) -> list[str]:
    out = []
    for k, t in res["tasks"].items():
        line = (f"    {k}: paired mean {t['mean_gain']:+.2f} (seed SD of the paired difference {t['sd']:.2f}), ahead "
                f"in {t['seeds_ahead']} of {len(SEEDS)} seeds, largest single loss {t['largest_loss']:.2f}, identical "
                f"on {t['identical_seeds']} of {len(SEEDS)} seeds")
        if t["idle"]:
            line += (f"; herding idle on {k} (the draw repeats rows or leaves nothing to swap in; the herded twin is "
                     "its base bit for bit): herding is not read there")
        out.append(line)
    return out


def _verdict(p: dict) -> str:
    if p["outcome"] == "not readable":
        return (f"is not readable under the adoption rule: herding is idle on {', '.join(p['idle'])}, and the rule "
                f"needs a pass on {p['need']} of {len(p['tasks'])} counted backtests (keep_change, deciding nothing "
                f"here: {p['reason']})")
    return f"{p['outcome']} the adoption rule ({p['reason']})"


def text(res: dict) -> str:
    if res["void"]:
        return f"{res['board']} ({res['test_board']}): VOID: " + "; ".join(res["why"]) + "\n"
    e, h = res["entry"], res["herded"]
    st = res["statement"]
    lines = [f"Reading of herded draws on the {res['board']}'s first-upload entry ({res['test_board']}; run "
             f"{res['runs']}, commit {res['commit'][:7]}): {h} against {e} on {len(res['counted'])} counted "
             f"backtests ({', '.join(res['counted'])}); the rule needs a pass on {st['need']} of "
             f"{len(res['counted'])}",
             f"  the herded entry {_verdict(st)}", *_rows(st)]
    for k, by_arm in res["reproduction"].items():
        for a, r in by_arm.items():
            if r["checked"]:
                lines.append(f"  {a} on {k}: reproduces its committed scores ({r['file']}, commit "
                             f"{str(r['commit'])[:7]}), largest difference {r['max_abs']:.1e} (void above {TOL:g})")
            elif a == e:
                lines.append(f"  {a} on {k}: the reproduction could not be checked (no committed harness result holds "
                             "its scores)")
            else:
                lines.append(f"  {a} on {k}: no committed harness result holds its scores (not checked)")
    b, c = res["pair"]
    if res["pair_is_the_statement"]:
        lines.append(f"  the registered pair ({c} against {b}) is the entry and its herded twin themselves: reported "
                     "above, it decides nothing more")
    else:
        rp = res["registered_pair"]
        word = {"passes": "would be kept", "does not pass": "would not be kept",
                "not readable": "is not readable"}[rp["outcome"]]
        lines.append(f"  the registered pair, deciding nothing: {c} against {b} {word} "
                     + (f"(herding idle on {', '.join(rp['idle'])}; {rp['reason']})" if rp["idle"] else
                        f"({rp['reason']})"))
        lines += _rows(rp)
    cols = (*POSITION_ONLY, SHAPE_GROUP)
    lines.append("  position-only metrics, descriptive and deciding nothing (" + ", ".join(POSITION_ONLY)
                 + f"; {SHAPE_GROUP} skill): mean over seeds [seed-to-seed SD]; for a herded arm, its per-seed "
                 "difference from the entry, mean [seed SD]. A herded draw is a different set of cells, so these "
                 "move with the change of draw alone:")
    for k, by_arm in res["position_only"].items():
        for a, v in by_arm.items():
            d = "" if a == e else " (minus the entry's: " + ", ".join(
                f"{v[m]['diff_mean']:+.4f} [{v[m]['diff_sd']:.4f}]" for m in cols) + ")"
            lines.append(f"    {k} {a}: " + ", ".join(f"{v[m]['mean']:.4f} [{v[m]['sd']:.4f}]" for m in cols) + d)
    lines.append("  " + LIMITS)
    if BOARDS[res["board"]].limits:
        lines.append("  " + BOARDS[res["board"]].limits)
    return "\n".join(lines) + "\n"


def _rel(p: Path) -> str:
    p = Path(p).resolve()
    root = HERE.parents[1]
    return str(p.relative_to(root)) if p.is_relative_to(root) else str(p)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("read")
    r.add_argument("--board", choices=sorted(BOARDS), required=True)
    r.add_argument("--runs", required=True, help="the harness run folder (task files) of the board")
    r.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    files = [out / f"n53_{a.board}.json", out / f"n53_{a.board}.txt"]
    if any(f.exists() for f in files):
        raise SystemExit(f"{out} already holds this board's reading: refusing to overwrite")
    res = reading(a.board, Path(a.runs))
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float) + "\n")
    files[1].write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
