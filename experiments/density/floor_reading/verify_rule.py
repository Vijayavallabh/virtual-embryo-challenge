# Since decision log 2026-10-21 (c) n_cells_from_decision applies the next count up (C); this checks rule 7 and G.
"""Cross-check the Monte Carlo's vectorised rule 7 and G reading against the repo's own code (read-only import).

For random 10-seed score tables it compares: the winner (vec.t2_build.count_override) and, in the floor case, the
level the code's ranking walks to (vec.t2_build._passing_ranking, which n_cells_from_decision uses: first level in
rule 7's order that the board allows, i.e. reading G), with the simulation's own vectorised implementation.
Run from the repo root: uv run --no-sync python -I <this file>
"""

import sys

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import sim_floor as S  # noqa: E402  (scratch module; its main() is not run on import)

from vec.t2_build import _passing_ranking, count_override  # noqa: E402

rng = np.random.default_rng(7)
n = 4000
true_gain = rng.normal(0.6, 0.6, (n, 3))
Z = rng.standard_normal((n, 10, 4)) @ np.linalg.cholesky(S.noise_cov("fit", 0.94)).T
d = true_gain[:, None, :] + Z[:, :, :3] - Z[:, :, 3:4]
base = 90 + rng.normal(0, 1, (n, 10))
labels = ["0.25max", "0.5max", "0.75max"]
mism_w = mism_g = n_floor = n_dis = 0
for i in range(n):
    scores = {"1max": {s: float(base[i, s]) for s in range(10)}}
    for j, lv in enumerate(labels):
        scores[lv] = {s: float(base[i, s] + d[i, s, j]) for s in range(10)}
    override, rec = count_override(scores)
    # simulation's vectorised rule on the same differences (recomputed from the rounded floats the code saw)
    dd = np.array([[scores[lv][s] - scores["1max"][s] for lv in labels] for s in range(10)])
    gain = dd.mean(0)
    ahead = (dd > 0).sum(0)
    passing = (gain > 0.3) & (ahead >= 8)
    w = np.where(passing, gain, -np.inf)
    win = (2 - int(np.argmax(w[::-1]))) if passing.any() else 3
    code_win = 3 if override is None else labels.index(override)
    mism_w += win != code_win
    if code_win == 0:
        n_floor += 1
        ranking = _passing_ranking({"record": rec}, override, "heart")
        allowed = [lv for lv in ranking if lv != "0.25max"]  # 0.25max is the level below the floor
        code_G = labels.index(allowed[0]) if allowed else 3
        both = passing[1] and passing[2]
        sim_G = (2 if gain[2] >= gain[1] else 1) if both else (1 if passing[1] else (2 if passing[2] else 3))
        mism_g += code_G != sim_G
        n_dis += both and gain[2] >= gain[1]
print(f"tables {n}: winner mismatches {mism_w}; floor cases {n_floor}; G-reading mismatches {mism_g}; "
      f"disagreement cases {n_dis}")
