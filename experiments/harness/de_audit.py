#!/usr/bin/env python
"""de_score audit of harness results (rows written since ab072da keep each candidate's metric panel).

Per pseudo-task and candidate, means over seeds of de_score's parts (veckit, pinned git cc3859a, common/core_metrics):
  raw      the share of the truth's DE genes (|change| >= 0.25 against the reference's 10% working copy) that the
           prediction's own ranking of its change puts in the matching slot (its top n_up, its bottom n_dn)
  chance   the same overlap for a ranking by the reference's own expression, the better of its two signs (the
           "expression-magnitude null")
  score    (raw - chance) / (1 - chance); a prediction with no change (the floor) scores exactly 0
and de_direction, the DE-set sizes, and each DE metric's board points against the floor, weight x (skill - 0.5) x 100
(T1 weights 0.25 and 0.25; T2 0.125 each). Ahead = seeds where raw > chance.

Usage: uv run python experiments/harness/de_audit.py RESULT.json [...] > audit.txt
"""

from __future__ import annotations

import json
import sys

import numpy as np

from vec.boards import GROUPS


def weights(task: str) -> dict[str, float]:
    return {m: g * w for _, g, ms in GROUPS[task] for m, w in ms}


def audit(path: str) -> list[str]:
    r = json.load(open(path))
    task, rows = r["task"]["task"], r["rows"]
    if not all("metrics" in row for row in rows):
        return [f"{r['task']['key']}: rows without metric panels (written before ab072da)"]
    w = weights(task)
    out = [f"{r['task']['key']} ({len(rows)} seeds, commit {r['provenance']['commit'][:7]}, truth "
           f"{r['board']['truth_cells']} cells, DE reference {r['board']['ref_cells']} cells, candidates "
           f"{rows[0]['n_cells']})",
           f"  {'candidate':18s} {'raw':>6s} {'chance':>6s} {'score':>7s} {'raw>chance':>10s} {'n_up':>5s} "
           f"{'n_dn':>5s} {'de_dir':>7s} {'pts de_score':>12s} {'pts de_dir':>10s} {'board':>6s}"]
    for c in r["candidates"]:
        m = [row["metrics"][c] for row in rows]
        raw = np.array([x["_de_raw"] if x["_de_raw"] is not None else 0.0 for x in m], dtype=float)
        ch = np.array([x["_de_chance"] for x in m], dtype=float)
        sc = np.array([x["de_score"] for x in m], dtype=float)
        dd = np.array([x["de_direction"] for x in m], dtype=float)
        sk = {k: np.mean([row["skills"][c][k] for row in rows]) for k in ("de_score", "de_direction")}
        pts = {k: w[k] * (sk[k] - 0.5) * 100 for k in sk}
        raw_txt = "  none" if c == "floor" else f"{raw.mean():6.3f}"
        ahead = "-" if c == "floor" else f"{int((raw > ch).sum())}/{len(rows)}"
        out.append(f"  {c:18s} {raw_txt} {ch.mean():6.3f} {sc.mean():+7.3f} {ahead:>10s} "
                   f"{np.mean([x['_n_up'] for x in m]):5.0f} {np.mean([x['_n_dn'] for x in m]):5.0f} "
                   f"{dd.mean():+7.3f} {pts['de_score']:+12.2f} {pts['de_direction']:+10.2f} "
                   f"{np.mean([row['scores'][c] for row in rows]):6.2f}")
    return out


def main(argv=None) -> None:
    paths = (argv if argv is not None else sys.argv[1:])
    for p in paths:
        print("\n".join(audit(p)))
        print()


if __name__ == "__main__":
    main()
