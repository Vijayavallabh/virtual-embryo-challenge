#!/usr/bin/env python3
"""The forecast for the files as built, as a function of its judgement weights (2026-11-14).

Aggregation of committed practice results only, not a try: no target, no data file, no scorer, no build, no new score.
Every figure comes through experiments/gate/forecast.py's loaders and world table and forecast_built.py's file set
(rows R6, R7 and H11); nothing here adds a world or a score.

Why. On 2026-11-14 the captain asked that the organisers' logical replies be assumed for the questions still open
(docs/ORGANIZER_QUESTIONS.md, "Working answers (2026-11-14)"; docs/PLAN.md decision log 2026-11-14 (a)). Three of
the forecast's judgement weights lie near three of those replies, and none is defined by one:
  emb_p_deep     the chance that E7.75 shares its deeper bracket's expression profile, not that it shared a batch;
  T1 depth       the chance that E12.5 is as deep as E10.5's library 3 or deeper, not that its protocol was the same;
  t3_standin_A   the weight of a fold type (a Gata4 factor read on the other Gata4 replicate), not the chance that
                 the test knockout shared a slide with those replicates.
An independent check of 2026-11-14 said so in these words: the restated totals "should be presented as a sensitivity
of the A weight, not as what that reply implies". This table is that sensitivity. It answers one question only: how
far could any set of replies move the forecast, if each weight went to an end of its range.

Usage, from the repository root:
  uv run --frozen --no-sync python -I -B experiments/gate/forecast_weights.py \
      --json experiments/gate/results/forecast_weights_2026-11-14.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import forecast as F  # noqa: E402
import forecast_built as FB  # noqa: E402
import forecast_t1_pair as FP  # noqa: E402

E = F.E
LINES = (203.00, 205.32, 209.00, 215.41)  # the lines of forecast_built.py (screened read of 2026-11-13)


def depth_weights(p_deep: float, deep: dict[float, float] | None = None) -> dict[float, float]:
    """forecast.py's depth weights: the shallow side by its census split, the deep side by `deep` (default the
    census split of forecast.T30_DEEP_SPLIT)."""
    deep = deep or F.T30_DEEP_SPLIT
    return {**{d: (1 - p_deep) * f for d, f in F.T30_SHALLOW_SPLIT.items()}, **{d: p_deep * f for d, f in deep.items()}}


C = F.central()
# name: (the weights moved, the T1 depth weights or None for the committed split)
CASES = {
    "as committed": ({}, None),
    "embryo, emb_p_deep at its low 0.30": ({"emb_p_deep": 0.30}, None),
    "embryo, emb_p_deep at its high 0.60": ({"emb_p_deep": 0.60}, None),
    "embryo, emb_p_deep 0 (never its deeper bracket's profile)": ({"emb_p_deep": 0.0}, None),
    "embryo, emb_p_deep 1 (always its deeper bracket's profile)": ({"emb_p_deep": 1.0}, None),
    "T1, depth exactly library 3's": ({}, depth_weights(1.0, {1.0: 1.0, 1.37: 0.0})),
    "T1, depth 1.37 times library 3's": ({}, depth_weights(1.0, {1.0: 0.0, 1.37: 1.0})),
    "T1, shallower than library 3 for certain": ({}, depth_weights(0.0)),
    "T3, stand-in A at its high 0.75": ({"t3_standin_A": 0.75}, None),
    "T3, stand-in A at its low 0.25": ({"t3_standin_A": 0.25}, None),
    "T3, stand-in A 0.15": ({"t3_standin_A": 0.15}, None),
    "T3, stand-in A 0": ({"t3_standin_A": 0.0}, None),
    "T3, stand-in A 1": ({"t3_standin_A": 1.0}, None),
    "T3, stand-in A 0 and the factor's carry at its low 0.35": ({"t3_standin_A": 0.0, "t3_q_carry": 0.35}, None),
    "every weight at the end that suits the files least": (
        {"emb_p_deep": 0.0, "t3_standin_A": 0.0, "t3_q_carry": 0.0}, depth_weights(0.0)),
    "every weight at the end that suits the files most": (
        {"emb_p_deep": 1.0, "t3_standin_A": 1.0, "t3_q_carry": 1.0}, depth_weights(1.0, {1.0: 0.0, 1.37: 1.0})),
}


def table() -> dict:
    inp = F.Inputs()
    arms = FP.load_arms(inp, ("s1_base", "s1_cm60", *FB.BUILT_PAIR))
    shift = float(np.mean([np.mean(v) for v in (
        np.array([r["scores"]["b0_trend"] - r["scores"]["b0"] for r in json.loads(
            (F.REPO / FB.TREND.format(t)).read_text())["rows"]]) for t in FB.TREND_TASKS)]))
    out = {}
    for name, (moved, dw) in CASES.items():
        par = {**C, **moved}
        bd = F.boards(inp, "new", par)
        t1 = FP.board_t1(inp, par["t1_prior"], par["t1_p_deep"], par["t1_k"], arms, FB.BUILT_PAIR, depth_w=dw)
        heart = dict(bd["heart"])
        heart["a"], heart["b"] = heart["a"] + shift, heart["b"] + shift
        b = {**bd, "T1": t1, "heart": heart}
        tot, _ = F.draw(b)
        out[name] = {"weights_moved": moved, "t1_depth_weights": dw,
                     "mean": float(tot.mean()), "median": float(np.percentile(tot, 50)),
                     "p05": float(np.percentile(tot, 5)), "p95": float(np.percentile(tot, 95)),
                     "boards_better_of_two": {k: float(E.reading(b[k], "a")) for k in ("T1", "heart", "embryo", "T3")},
                     "p_at_least": {f"{c:.2f}": float(np.mean(tot >= c)) for c in LINES}}
    return {"label": "the forecast for the files as built as a function of its judgement weights; a sensitivity, "
                     "not what any organiser reply implies", "central_weights": C, "heart_shift": shift,
            "draws": int(F.N_DRAWS), "rng_seed": int(F.SEED), "cases": out}


def text(res: dict) -> str:
    lines = [res["label"], f"  {'case':<58}   mean  5th-95th      T1    T3  heart embryo  "
             + " ".join(f"P(>={c:g})" for c in LINES)]
    for name, v in res["cases"].items():
        b = v["boards_better_of_two"]
        lines.append(f"  {name:<58} {v['mean']:6.1f}  {v['p05']:5.1f}-{v['p95']:5.1f}  {b['T1']:5.1f} {b['T3']:5.1f}  "
                     f"{b['heart']:5.1f}  {b['embryo']:5.1f}  "
                     + " ".join(f"{100 * v['p_at_least'][f'{c:.2f}']:8.2f}%" for c in LINES))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", metavar="PATH")
    a = ap.parse_args(argv)
    res = table()
    print(text(res), end="")
    if a.json:
        p = Path(a.json)
        if p.exists():
            raise SystemExit(f"{p} exists: refusing to overwrite")
        p.write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
