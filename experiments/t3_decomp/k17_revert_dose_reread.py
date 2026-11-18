"""Try K17 (experiments/harness/TRIES.md, registered before any read): the exponent of T3's wild-type-only factor,
read from the two runs try K16 made, under a guard that tolerates what the machines do not reproduce.

Try K16 (k16_revert_dose.py) was void on both machines by its own guard: the base arm had to equal the committed
scores of tries K8b and K7b seed for seed to 1e-9. On the local machine, where K7 and K8 ran, the base equals K8b
exactly on g1 and g2 and differs on m in two seeds by at most 0.004 board points; on the DGX, where T3 builds draw a
few cells differently (docs/COMPUTE.md, "What reproduces where"), it differs in more seeds. Nothing of K16 was read
but those base-arm differences. This try reads the same task files; it builds nothing and scores nothing.

Arms, tasks, forms and seeds: K16's, unchanged (k16_revert_dose.ARMS, TASKS, SEEDS). The reading of one run is K16's
own `summarise` and `decide`, with one change to what makes it void:

  GUARD (each run)    in place of the 1e-9 reproduction: on every task, grp_e050's mean over the ten seeds is within
                      GUARD_MEAN (0.1 board points, a third of the reading's 0.3) of the committed base's mean (K8b on
                      g1 and g2, K7b on m), and no seed's base score is more than GUARD_SEED (0.5) from the committed
                      one. On the local run the base's distances were known before these limits were set, so there the
                      guard states a tolerance and tests nothing; on the DGX run they were not, so there it can void.
  FACTOR RECORD       each task file's factor record (the fitted factor's summary: no score) equals the one try K7
                      committed for that task: a score-free check of drift in code or data.
  The other void conditions are K16's (ten seeds, every arm, a clean tree, one commit per run), and the two runs must be
  at one commit.

The outcome, fixed in the row:
  REQUEST   exponent e, only if K16's `decide`, applied unchanged to each run, names the same e with its request true
            (carries, steps, factor-specific) on both. A priced question to the captain after the cold check row K16
            words.
  Else      no request; the second upload stays at exponent 0.5. That includes: either run void; either run naming no
            exponent or failing the factor-specific check; the two runs naming different exponents.

  read --local DIR --dgx DIR --out DIR     apply the reading; write reading.json and reading.txt
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k16_revert_dose as k16  # noqa: E402  (the arms, the tasks, summarise, decide, the text of one run's reading)

GUARD_MEAN = 0.1  # board points: grp_e050's ten-seed mean against the committed base's mean, per task
GUARD_SEED = 0.5  # board points: any one seed's base score against the committed one


def base_distance(name: str, rows: list[dict], refs: Path) -> dict | None:
    """grp_e050 against the committed base of the task (K8b on g1 and g2, K7b on m): the difference of the ten-seed
    means, the largest absolute per-seed difference and the seeds that differ by more than 1e-9. None when the
    committed file or one of its seeds is missing."""
    if not rows or any(k16.BASE not in q["scores"] for q in rows):
        return None
    folder, ref_arm = k16.REPRO[name]
    p = refs / folder / "tasks" / f"{name}.json"
    if not p.exists():
        return None
    ref = {x["seed"]: x["scores"].get(ref_arm) for x in json.loads(p.read_text())["rows"]}
    if any(ref.get(q["seed"]) is None for q in rows):
        return None
    d = np.array([q["scores"][k16.BASE] - ref[q["seed"]] for q in rows], dtype=float)
    return {"mean_diff": float(d.mean()), "max_abs_diff": float(np.abs(d).max()),
            "seeds_off": [q["seed"] for q, v in zip(rows, d, strict=True) if abs(v) > k16.REPRO_TOL]}


def guard(name: str, rows: list[dict], refs: Path) -> str | None:
    """None when the base passes the guard; else why not."""
    d = base_distance(name, rows, refs)
    if d is None:
        return f"{name}: the committed base is missing or lacks a seed"
    if abs(d["mean_diff"]) > GUARD_MEAN or d["max_abs_diff"] > GUARD_SEED:
        return (f"{name}: {k16.BASE} is off the committed base by {d['mean_diff']:+.3f} in the mean and "
                f"{d['max_abs_diff']:.3f} in one seed (limits {GUARD_MEAN} and {GUARD_SEED})")
    return None


FACTOR_REF = "k7_wt_revert"  # the committed factor record of each task (try K7's task files)


def factor_drift(name: str, task: dict, refs: Path) -> str | None:
    """None when the task file's factor record equals the one try K7 committed for the task; else why not. The record
    summarises the fitted factor (its groups, spreads, clipped genes and samples) and holds no score."""
    p = refs / FACTOR_REF / "tasks" / f"{name}.json"
    if not p.exists():
        return f"{name}: the committed factor record is missing"
    want = json.loads(p.read_text()).get("factor_record")
    if want is None or task.get("factor_record") != want:
        return f"{name}: the factor record differs from the one try K7 committed"
    return None


def read_run(runs: Path, refs: Path = k16.REF_ROOT) -> dict:
    """K16's reading of one run folder with its 1e-9 reproduction replaced by `guard`, plus the factor-record check.
    Every other void condition, `summarise` and `decide` are K16's own code. Also returned: the base's distance from
    the committed base per task, and the sha256 of each task file read."""
    saved = k16.reproduction
    k16.reproduction = guard
    try:
        res = k16.reading(runs, refs)
    finally:
        k16.reproduction = saved
    dist, sums, drift = {}, {}, []
    for name in k16.TASKS:
        p = runs / f"{name}.json"
        if not p.exists():
            continue
        raw = p.read_bytes()
        sums[name] = hashlib.sha256(raw).hexdigest()
        task = json.loads(raw)
        dist[name] = base_distance(name, sorted(task.get("rows", []), key=lambda x: x["seed"]), refs)
        miss = factor_drift(name, task, refs)
        if miss:
            drift.append(miss)
    if drift:
        res = {"void": True, "why": [*res["why"], *drift]}
    return res | {"base_distance": dist, "task_sha256": sums}


def combine(local: dict, dgx: dict) -> dict:
    """The outcome fixed in row K17."""
    if local["void"]:
        return {"request": False, "exponent": None, "why": "the local run is void: " + "; ".join(local["why"])}
    e = local["named_exponent"]
    if e is None or not local["request"]:
        why = "the local run names no exponent" if e is None else "the local run fails the factor-specific check"
        return {"request": False, "exponent": None, "why": why}
    if dgx["void"]:
        return {"request": False, "exponent": None, "why": "the DGX run is void: " + "; ".join(dgx["why"])}
    if local["commit"] != dgx["commit"]:
        return {"request": False, "exponent": None,
                "why": f"the two runs are at different commits ({local['commit']}, {dgx['commit']})"}
    d = dgx["named_exponent"]
    if d is None or not dgx["request"]:
        why = "the DGX run names no exponent" if d is None else "the DGX run fails the factor-specific check"
        return {"request": False, "exponent": None, "why": why}
    if d != e:
        return {"request": False, "exponent": None,
                "why": f"the runs name different exponents (local {e:g}, DGX {d:g})"}
    return {"request": True, "exponent": e, "why": f"both runs name exponent {e:g} with the request true"}


def reading(local_runs: Path, dgx_runs: Path, refs: Path = k16.REF_ROOT) -> dict:
    local = read_run(local_runs, refs)
    dgx = read_run(dgx_runs, refs)
    out = combine(local, dgx)
    out["outcome"] = (f"exponent {out['exponent']:g}: a priced question to the captain after the cold check row K16 "
                      "words (T3's frozen second-upload rule; a builder option and a rehearsal build); the try changes "
                      "no entry" if out["request"] else "no request; the second upload stays at exponent 0.5")
    return {"local": local, "dgx": dgx} | out


def _dist(d: dict) -> str:
    def one(v) -> str:
        if v is None:
            return "n/a"
        return (f"mean {v['mean_diff']:+.4f}, largest seed {v['max_abs_diff']:.4f}, "
                f"{len(v['seeds_off'])} of {len(k16.SEEDS)} seeds differ")

    return "; ".join(f"{n} {one(v)}" for n, v in d.items())


def text(res: dict) -> str:
    lines = [f"Try K17: {res['outcome']} ({res['why']}).",
             "the base against the committed base, local run: " + _dist(res["local"]["base_distance"]),
             "the base against the committed base, DGX run: " + _dist(res["dgx"]["base_distance"]),
             "task files read (sha256): " + "; ".join(
                 f"{m} {n} {h[:16]}" for m in ("local", "dgx") for n, h in res[m]["task_sha256"].items()),
             "", "LOCAL RUN (K16's reading under the guard):",
             k16.text(res["local"]).rstrip("\n").replace("Try K16", "run"),
             "", "DGX RUN (K16's reading under the guard):",
             k16.text(res["dgx"]).rstrip("\n").replace("Try K16", "run")]
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    res = reading(Path(a.local), Path(a.dgx))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("read")
    r.add_argument("--local", required=True)
    r.add_argument("--dgx", required=True)
    r.add_argument("--out", required=True)
    r.set_defaults(func=cmd_read)
    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
