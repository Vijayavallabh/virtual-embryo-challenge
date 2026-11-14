#!/usr/bin/env python
"""Reading of try E43 of experiments/harness/TRIES.md (registered before any run): `fgw_partial_grp`, the partial
coupling of try E41b with its transported share taken from the lineage-group labels in place of the cell-type
labels, scored on the nine practice tasks by `n47_embryo_h_arms.py` on the branch `campaign/emb-h`.

Why. Row E41 named, before its run, one way a pass of `fgw_partial` could be a fit and not an effect: "a gain that
depends on the label level". Cell-type vocabularies are not harmonised across stages and `docs/DATA_NOTES.md`
prescribes the merged lineage level for cross-stage work. A cold adversarial check of the pass asked for this read.

What it reads. Only per-seed board scores committed by counted tries: the stored runs of E24, E25 and E28, the
committed task files of tries E41 (`results/n45_embryo_f_arms/tasks`: `fgw_partial`) and the committed task files
of try E43 given by `--runs`. It opens no data file, runs no candidate and no scorer.

Validity, a check that can fail. The run rescored `fgw`, `fgw_partial` and the floor; each must agree with its
committed scores seed for seed within nine_tasks.TOL on every task, or the reading is void.

Reading, fixed in the row before the run: the adoption rule (vec.adoption.keep_change) for `fgw_partial_grp`
against `fgw` on the two counted embryo backtests, and row N33's first-upload rule beside the second upload as
built. The partial form is "label-robust" if the arm passes both, as `fgw_partial` did. If it does not, row E41's
clause applies: the gain depends on the label level, and the partial form is not put forward. What a label-robust
form leads to is fixed in row E43 together with row E42's reading; this script requests nothing.
Reported, deciding nothing: the arm's own mean minus `fgw`'s and minus `fgw_partial`'s on every task.

  python -I -B experiments/harness/n47_embryo_h_readings.py --runs DIR --out DIR
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n37_first_upload_new_arms as n37  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec.adoption import keep_change  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ENTRY = "fgw"  # the first upload as built
SECOND = nt.RA  # fgw_aligned_series, the second upload as built
PARTIAL = "fgw_partial"  # try E41b, share from the cell-type labels
ARM = "fgw_partial_grp"  # the same arm, share from the lineage-group labels
RESCORED = ("floor", ENTRY, PARTIAL)  # scored again by the run; must match the record
BACKTESTS = ("e75", "e725")  # the two counted embryo backtests among the nine tasks
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
N45_RUNS = HERE / "results" / "n45_embryo_f_arms" / "tasks"


def adoption(scores: dict, tasks: dict, arm: str, base: str = ENTRY) -> dict:
    """The adoption rule on the two counted backtests, seed for seed (rows are in seed order)."""
    paired = {}
    for n in BACKTESTS:
        s = scores[tasks[n].key]
        paired[n] = {i: (float(b), float(a)) for i, (b, a) in enumerate(zip(s[base], s[arm], strict=True))}
    return keep_change(paired)


def rescored(committed: dict, run: dict, tasks: dict) -> tuple[dict, list[str]]:
    """Per rescored arm, the largest absolute difference between the run's per-seed scores and the committed ones over
    the nine tasks, and the reasons the run is void (an arm missing from the run, or off by more than nt.TOL)."""
    worst, why = {}, []
    for arm in RESCORED:
        big = 0.0
        for n, t in tasks.items():
            new, old = run.get(t.key, {}).get(arm), committed.get(t.key, {}).get(arm)
            if new is None:
                why.append(f"{n}: the run did not rescore {arm}")
                continue
            if old is None:
                why.append(f"{n}: no committed scores for {arm}")
                continue
            d = float(np.abs(np.asarray(new, float) - np.asarray(old, float)).max())
            big = max(big, d)
            if not d <= nt.TOL:
                why.append(f"{n}: {arm} differs from its committed scores by {d:.2e}")
        worst[arm] = big
    return worst, why


def check(scores: dict, tasks: dict) -> dict:
    names = [*nt.DEEP, *nt.OTHER]
    readable = all(ARM in scores[tasks[n].key] and np.isfinite(scores[tasks[n].key][ARM]).all() for n in names)

    def minus(a: str, b: str) -> dict:
        return {n: float(scores[tasks[n].key][a].mean() - scores[tasks[n].key][b].mean()) for n in names}

    first = n37.price(scores, tasks, candidates=(ARM, PARTIAL))
    out = {"arm": ARM, "readable": readable,
           "own_minus_fgw": minus(ARM, ENTRY) if readable else None,
           "own_minus_partial": minus(ARM, PARTIAL) if readable else None,
           "adoption": adoption(scores, tasks, ARM) if readable else None,
           "adoption_partial": adoption(scores, tasks, PARTIAL),
           "first_upload": first}
    out["label_robust"] = bool(readable and out["adoption"]["keep"] and first["first_uploads"][ARM]["passes"])
    return out


def reading(run_dirs) -> dict:
    committed, why, meta = nt._load([*nt.STORED, E28_RUNS, N45_RUNS])
    run, why_run, meta_run = nt._load(list(run_dirs))
    why, meta = [*why, *why_run], sorted({*meta, *meta_run})
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (ENTRY, SECOND, PARTIAL):
            if arm not in committed.get(t.key, {}):
                why.append(f"{n}: no committed scores for {arm}")
        if t.key not in run:
            why.append(f"{n}: no task file in the run")
    worst, why_rescored = rescored(committed, run, tasks)
    why = [*why, *why_rescored]
    if why:
        return {"void": True, "why": sorted(set(why)), "rescored_max_abs_difference": worst}
    scores = {t.key: {**committed[t.key], **({ARM: run[t.key][ARM]} if ARM in run[t.key] else {})}
              for t in tasks.values()}
    return {"void": False, "why": [], "runs": meta, "rescored_max_abs_difference": worst, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    lines = [f"Reading of try E43 (runs {', '.join(res['runs'])}): the partial coupling with its share from the "
             "lineage-group labels, on the nine practice tasks; aggregation of committed scores",
             "the partial form is " + ("label-robust: the arm passes both rules" if res["label_robust"] else
                                       "NOT label-robust: the arm does not pass both rules (row E41's clause applies)"),
             "  rescored arms against their committed scores, largest difference over tasks and seeds: "
             + ", ".join(f"{a} {res['rescored_max_abs_difference'][a]:.1e}" for a in RESCORED)
             + f" (void above {nt.TOL})"]
    if not res["readable"]:
        return "\n".join([*lines, f"  {ARM}: not readable"]) + "\n"
    for arm, a in ((ARM, res["adoption"]), (PARTIAL, res["adoption_partial"])):
        parts = [f"{n} {t['mean_gain']:+.2f} ({t['seeds_ahead']}/10)" for n, t in a["tasks"].items()]
        v = res["first_upload"]["first_uploads"][arm]
        r = v["rule"]
        lines.append(f"  {arm:<16} adoption against {ENTRY}: " + ", ".join(parts)
                     + f" -> {'kept' if a['keep'] else 'not kept'}; first upload S(0.45) {v['S0.45']:+.2f}, clauses "
                     + " ".join("y" if r[k] else "n" for k in sorted(r)))
    lines.append(f"  {ARM} own mean minus {ENTRY}: " + " ".join(f"{n} {res['own_minus_fgw'][n]:+.2f}" for n in names))
    lines.append(f"  {ARM} own mean minus {PARTIAL}: "
                 + " ".join(f"{n} {res['own_minus_partial'][n]:+.2f}" for n in names))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of try E43")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n47.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n47.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
