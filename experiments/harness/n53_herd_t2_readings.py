#!/usr/bin/env python
"""Reading of tries H24 (heart) and E46 (embryo) of experiments/harness/TRIES.md: herded draws on each T2 board's
first-upload entry, scored by `experiments/harness/run.py run` on the board's counted backtests (decision log
2026-11-18 (c) fixes its use before any run).

Why. Try T5 read herded draws on Task 1 (within each category's quota a drawn cell is swapped for an undrawn one when
that brings the sample's mean expression closer to the pool's; `vec.core.herd_rows`). Rows H4 and E9 registered the
same option on the T2 boards (`b0_herd`, `fgw_herd`). This reading puts the option on each board's first-upload entry
as the harness builds it: the heart's `b0_trend` against `b0_trend_herd` (`b0_trend` with the draw herded as
`b0_herd`'s), the embryo's `fgw` against `fgw_herd` (the embryo's first upload is `fgw` itself: `build_test.py embryo
--method fgw --expression mixture`, so its herded twin is the registered `fgw_herd`).

What it reads. Only the harness task files of one board in `--runs` (per-seed scores, skills and metrics), and the
entry's committed task files (`COMMITTED`) to check that the run reproduces them. It opens no data file, runs no
candidate and no scorer.

Validity, checks that can fail (otherwise the reading is void and decides nothing). (1) The run's task files pass
vec.harness.check_results: one commit, a clean tree, the current config, the bake protocol, seeds 0-9 and their build
seeds; no backtest failed. (2) Every counted backtest of the board (configured, counted, and every arm applies) holds
the four arms (the entry, its herded twin and the registered pair) with a finite score on each of the ten seeds. (3)
Where a committed harness result holds the entry's scores on a backtest (`COMMITTED`), the run's entry agrees with them
within TOL on every seed; where none does, the reading says the reproduction could not be checked.

Statement, per board: vec.adoption.keep_change of the herded entry against the entry on the counted backtests ("the
herded entry passes the adoption rule" or "does not pass"), with each backtest's paired mean, seeds ahead and largest
single loss. The registered pair (b0 / b0_herd, fgw / fgw_herd) is reported the same way and decides nothing. Reported
too, deciding nothing: the position-only metrics (the shape group: d2_shape, occupancy_dice, scale_log_ratio) and the
shape group's skill of every arm beside the entry's (herding must not depend on position; a scored arm passed the T2
head's cloud guard, which raises otherwise).

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

TOL = 1e-9  # per seed, board points: the entry's reproduction of its committed scores
SHAPE_GROUP = "shape_growth"
POSITION_ONLY = tuple(m for g, _, ms in GROUPS["T2"] if g == SHAPE_GROUP for m, _ in ms)


@dataclass(frozen=True)
class Spec:
    test_board: str
    entry: str  # the board's first-upload entry as the harness builds it
    herded: str  # its herded twin
    pair: tuple[str, str]  # the registered pair (rows H4 and E9): base, herded
    committed: Path  # the committed harness task files that hold the entry's scores

    @property
    def arms(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys((self.entry, self.herded, *self.pair)))


BOARDS = {
    "heart": Spec("T2:heart:test_extrap", "b0_trend", "b0_trend_herd", ("b0", "b0_herd"),
                  HERE / "results" / "trend_heart_v1" / "tasks"),  # row H11's run of b0_trend
    "embryo": Spec("T2:embryo:test_interp", "fgw", "fgw_herd", ("fgw", "fgw_herd"),
                   HERE / "results" / "emb_fgw_v1" / "tasks"),  # the embryo entry's run of 2026-11-01
}

LIMITS = ("Limits: herding has no biological content (it moves the drawn sample's mean panel expression towards the "
          "pool's and nothing else); the scorer's own reference working copy carries much of the de terms' noise, "
          "which herding the prediction cannot remove; a pass on two or three backtests of one board is in-sample for "
          "nothing but is still few tasks.")


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
    """Per counted backtest, whether the run's entry reproduces its committed scores: {"checked": bool, "max_abs":
    float | None, "file": str | None}; and the reasons the run is void (a difference above TOL)."""
    committed = {}
    if Path(committed_dir).is_dir():
        for p in sorted(Path(committed_dir).glob("*.json")):
            r = json.loads(p.read_text())
            if isinstance(r, dict) and "task" in r and "rows" in r:
                committed[r["task"]["key"]] = (r, p)
    out, why = {}, []
    for k in keys:
        have = committed.get(k)
        old = scores_of(have[0], spec.entry) if have else None
        if not complete(old):
            out[k] = {"checked": False, "max_abs": None, "file": None}
            continue
        new = scores_of(run[k], spec.entry)
        d = max(abs(new[s] - old[s]) for s in SEEDS)
        out[k] = {"checked": True, "max_abs": d, "file": _rel(have[1])}
        if not d <= TOL:
            why.append(f"{k}: {spec.entry} differs from its committed scores ({_rel(have[1])}) by {d:.2e}")
    return out, why


def paired(run: dict[str, dict], keys: list[str], base: str, change: str) -> dict:
    """keep_change's verdict of `change` against `base`, with each backtest's largest single loss (the largest
    per-seed base - change, 0 if the change is ahead or level on every seed)."""
    pairs = {k: {s: (scores_of(run[k], base)[s], scores_of(run[k], change)[s]) for s in SEEDS} for k in keys}
    out = keep_change(pairs)
    for k in keys:
        out["tasks"][k]["largest_loss"] = max(0.0, max(b - c for b, c in pairs[k].values()))
    return {"base": base, "change": change, **out}


def position(run: dict[str, dict], keys: list[str], arms) -> dict:
    """Per backtest and arm, the mean over seeds of each position-only metric (raw values) and of the shape group's
    skill (the group's within-group weights on the skills)."""
    weights = dict(next(ms for g, _, ms in GROUPS["T2"] if g == SHAPE_GROUP))
    out = {}
    for k in keys:
        rows = run[k]["rows"]
        out[k] = {a: {**{m: float(np.mean([row["metrics"][a][m] for row in rows])) for m in POSITION_ONLY},
                      SHAPE_GROUP: float(np.mean([sum(w * row["skills"][a][m] for m, w in weights.items())
                                                  for row in rows]))}
                  for a in arms}
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
    return {"void": False, "why": [], **head, "commit": commit, "passes": bool(stmt["keep"]), "statement": stmt,
            "registered_pair": pair, "pair_is_the_statement": tuple(spec.pair) == (spec.entry, spec.herded),
            "position_only": position(run, keys, spec.arms)}


def _rows(res: dict) -> list[str]:
    return [f"    {k}: paired mean {t['mean_gain']:+.2f}, ahead in {t['seeds_ahead']} of {len(SEEDS)} seeds, largest "
            f"single loss {t['largest_loss']:.2f}" for k, t in res["tasks"].items()]


def text(res: dict) -> str:
    if res["void"]:
        return f"{res['board']} ({res['test_board']}): VOID: " + "; ".join(res["why"]) + "\n"
    e, h = res["entry"], res["herded"]
    lines = [f"Reading of herded draws on the {res['board']}'s first-upload entry ({res['test_board']}; run "
             f"{res['runs']}, commit {res['commit'][:7]}): {h} against {e} on {len(res['counted'])} counted "
             f"backtests ({', '.join(res['counted'])})",
             f"  the herded entry {'passes' if res['passes'] else 'does not pass'} the adoption rule "
             f"({res['statement']['reason']})", *_rows(res["statement"])]
    for k, r in res["reproduction"].items():
        lines.append(f"  {e} on {k}: " + (f"reproduces its committed scores ({r['file']}), largest difference "
                                          f"{r['max_abs']:.1e} (void above {TOL:g})" if r["checked"] else
                                          "the reproduction could not be checked (no committed harness result holds "
                                          "its scores)"))
    b, c = res["pair"]
    if res["pair_is_the_statement"]:
        lines.append(f"  the registered pair ({c} against {b}) is the entry and its herded twin themselves: reported "
                     "above, it decides nothing more")
    else:
        lines.append(f"  the registered pair, deciding nothing: {c} against {b} "
                     f"{'would be kept' if res['registered_pair']['keep'] else 'would not be kept'} "
                     f"({res['registered_pair']['reason']})")
        lines += _rows(res["registered_pair"])
    lines.append("  position-only metrics, mean over seeds (" + ", ".join(POSITION_ONLY) + f"; {SHAPE_GROUP} skill), "
                 "deciding nothing:")
    for k, by_arm in res["position_only"].items():
        for a, v in by_arm.items():
            d = "" if a == e else " (minus the entry's: " + ", ".join(
                f"{v[m] - by_arm[e][m]:+.4f}" for m in (*POSITION_ONLY, SHAPE_GROUP)) + ")"
            lines.append(f"    {k} {a}: " + ", ".join(f"{v[m]:.4f}" for m in (*POSITION_ONLY, SHAPE_GROUP)) + d)
    lines.append("  " + LIMITS)
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
