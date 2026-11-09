#!/usr/bin/env python
"""Row E14's smoke abort (experiments/harness/TRIES.md; T2 embryo, candidate fgw_latent, vec.t2_latent), applied to the
harness run of seeds 0 and 1 on the E7.5 backtest (T2e/E7.25+E8.0->E7.5) before the full run.

The rule, fixed in row E14 before any run: the try stops if fgw_latent's mmd_u skill or its variogram skill is, on
average over seeds 0 and 1, more than 0.05 below fgw's (MAX_MEAN_DROP). An aborted run still counts as the try; a
run that continues is judged by the adoption rule on the full 10 seeds of both backtests (`run.py decide`).

It reads only the smoke run's own JSON (`run.py run --task T2e_E7.25_E8.0_E7.5 --candidates floor,fgw,fgw_latent
--seeds 0-1`): the two candidates' skills on two metrics, the run's commit and its clean flag. It refuses a run from a
dirty tree or with other seeds or candidates. It prints and writes only skills, their differences and the verdict: no
cell count, no test-board figure, no anchor. It refuses to overwrite its outputs (a rerun needs a new TRIES row).

Outputs, in experiments/harness/results/e14_latent_smoke/: e14_smoke.json and e14_smoke.txt, each stamped with
vec.provenance (commit, clean-tree flag, UTC, the run file's sha256).

Exit status: 0 = continue to the full run, 3 = abort.

Usage (DGX, repo root; any of the project's environments, it needs neither torch nor the data):
  uv run --frozen --no-sync python -I -B experiments/harness/e14_smoke.py --run runs/harness/e14_smoke
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import textwrap
import time
from pathlib import Path

from vec import harness, provenance

OUT = Path(__file__).resolve().parent / "results" / "e14_latent_smoke"
TASK_KEY = "T2e/E7.25+E8.0->E7.5"
RUN_FILE = "T2e_E7.25_E8.0_E7.5.json"  # run.py's safe(key) of TASK_KEY
BASE, CANDIDATE = "fgw", "fgw_latent"
METRICS = ("mmd_u", "variogram")
SEEDS = (0, 1)
MAX_MEAN_DROP = 0.05  # row E14: the mean skill difference (candidate - fgw) may not fall below -this
ABORT_EXIT = 3
NAMES = {"json": "e14_smoke.json", "text": "e14_smoke.txt"}  # write order
LABEL = "E14's smoke abort (fgw_latent against fgw, seeds 0 and 1 of the E7.5 backtest)"


def smoke_rule(rows: list[dict]) -> dict:
    """E14's smoke rule on the run's per-seed rows (each with `seed` and `skills`, {candidate: {metric: skill}}):
    for mmd_u and variogram the mean over seeds 0 and 1 of (fgw_latent's skill - fgw's); it aborts if either mean is
    below -MAX_MEAN_DROP. A non-finite skill makes the run invalid (ValueError), never a silent verdict."""
    by_seed = {int(r["seed"]): r for r in rows}
    if sorted(by_seed) != list(SEEDS) or len(rows) != len(SEEDS):
        raise ValueError(f"the smoke run must hold seeds {list(SEEDS)} once each, not {sorted(by_seed)}")
    checks = {}
    for m in METRICS:
        per_seed = []
        for s in SEEDS:
            skills = by_seed[s]["skills"]
            for name in (BASE, CANDIDATE):
                if name not in skills:
                    raise ValueError(f"seed {s}: the run has no {name}")
            a, b = float(skills[CANDIDATE][m]), float(skills[BASE][m])
            if not (math.isfinite(a) and math.isfinite(b)):
                raise ValueError(f"seed {s}: non-finite {m} skill")
            per_seed.append(a - b)
        mean = sum(per_seed) / len(per_seed)
        checks[m] = {"difference_per_seed": per_seed, "mean_difference": mean, "abort": bool(mean < -MAX_MEAN_DROP)}
    return {"checks": checks, "abort": any(c["abort"] for c in checks.values())}


def report_text(res: dict) -> str:
    verdict = "abort" if res["abort"] else "continue"
    lines = [f"{LABEL}: {verdict.upper()}"]
    lines += textwrap.wrap(
        f"Rule (row E14): the try stops if {CANDIDATE}'s mmd_u or variogram skill is, on average over seeds 0 and 1 of "
        f"the E7.5 backtest, more than {MAX_MEAN_DROP} below {BASE}'s. An aborted run still counts as the try.", 118)
    lines.append("")
    for m, ck in res["checks"].items():
        lines.append(f"  {m:<11s}skill difference ({CANDIDATE} - {BASE}) per seed "
                     + ", ".join(f"{d:+.4f}" for d in ck["difference_per_seed"])
                     + f"; mean {ck['mean_difference']:+.4f}: {'ABORT' if ck['abort'] else 'ok'}")
    s = res["provenance"]
    lines += ["", f"Run {res['run']['file']} sha256 {res['run']['sha256']}, commit {res['run']['commit']}, tree clean "
                  f"{res['run']['tree_clean']}", f"This check: commit {s['commit']}, tree clean {s['tree_clean']}, "
              f"UTC {s['utc']}", f"verdict: {verdict}"]
    return "\n".join(lines) + "\n"


def tree_state(exempt: list[Path]) -> dict:
    """HEAD and whether the tree is clean by vec.provenance's rule, the outputs exempt (None when git cannot tell)."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    return {"commit": commit, "clean": None if unclean is None else not unclean}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, help="the smoke run's folder (run.py run ... --seeds 0-1 --out)")
    ap.add_argument("--out", default=str(OUT), help="output folder (default: the registered one)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("E14's smoke verdict is already written; a rerun needs a new TRIES row")
    run_file = Path(a.run) / RUN_FILE
    res_run = json.loads(run_file.read_text())
    if res_run["task"]["key"] != TASK_KEY:
        sys.exit(f"{run_file}: not the {TASK_KEY} run")
    if res_run["seeds"] != list(SEEDS):
        sys.exit(f"{run_file}: the smoke run must hold seeds {list(SEEDS)} only")
    if not res_run["provenance"]["tree_clean"]:
        sys.exit(f"{run_file}: a run from a dirty tree does not decide the abort")
    missing = [c for c in (BASE, CANDIDATE) if c not in res_run["candidates"]]
    if missing:
        sys.exit(f"{run_file}: the run has no {missing}")
    paths = {k: out / v for k, v in NAMES.items()}
    exempt = [*paths.values(), *(provenance.stamp_path(p) for p in paths.values())]
    rule = smoke_rule(res_run["rows"])
    state = tree_state(exempt)
    res = {"label": LABEL, "task": TASK_KEY, "seeds": list(SEEDS), "max_mean_drop": MAX_MEAN_DROP,
           "checks": rule["checks"], "abort": rule["abort"],
           "run": {"file": RUN_FILE, "sha256": harness._sha256(run_file), "commit": res_run["provenance"]["commit"],
                   "tree_clean": res_run["provenance"]["tree_clean"]},
           "provenance": {"commit": state["commit"], "tree_clean": state["clean"],
                          "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    paths["json"].write_text(json.dumps(res, indent=1) + "\n")
    text = report_text(res)
    paths["text"].write_text(text)
    for p in paths.values():
        provenance.write_stamp(p)
    print(text, end="")
    return ABORT_EXIT if rule["abort"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
