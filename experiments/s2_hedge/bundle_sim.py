#!/usr/bin/env python
"""What the s2 set gains if final marks are one per set instead of per board (docs/S2_RULES.md, "Per-set marks").

Exact expectation over the committed per-seed payoffs, with no sampling:
- T1: the 60% s2 minus s1, over the T1 grid's 15 truths x 10 seeds.
- heart: s2 minus s1 over 3 transitions x 5 truths x 10 seeds.
- embryo: fgw_shape0 minus fgw over 2 pseudo-tasks x 10 seeds.

The boards are independent, and each T2 board counts half in the total (T1 + mean of the T2 boards). T3 is left out:
its s2 is a compliance hedge with a fixed order. "Own weight" w puts w on the released-like truth of each board (T1
bown_cown, heart own) and spreads the rest evenly over that board's other truths; "even" weights every truth alike.

Three sets are compared:
- A: the s2s in docs/S2_RULES.md (T1 at 60%, heart half trim, embryo fgw_shape0);
- B: T1 at 60%, heart H2 as first built, embryo redrawn (s1 at another seed, taken as a difference of 0);
- C: T1 at 60% only, heart and embryo redrawn.

Per-board marks: E[max(0, dT1)] + (E[max(0, dHeart)] + E[max(0, dEmbryo)]) / 2. Per-set marks:
E[max(0, dT1 + (dHeart + dEmbryo) / 2)].

Usage: uv run python experiments/s2_hedge/bundle_sim.py [--out experiments/s2_hedge/results/bundle_sim]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PER_SEED = ROOT / "experiments" / "s2_hedge" / "results" / "per_seed"
EMBRYO = ROOT / "experiments" / "harness" / "results" / "emb_shape_v1" / "tasks"
HEART_TAGS = ("heart-scope", "heart-scope-e825", "heart-scope-e875")


def _rows(pattern: str, root: Path = PER_SEED):
    for f in sorted(root.glob(pattern)):
        yield from json.loads(f.read_text())["rows"]


def deltas() -> dict:
    t1 = [(q, sc["bown_c0.6"] - sc["bown_cown"]) for r in _rows("t1-grid/*/t1-grid.json")
          for q, sc in r["scores"].items()]
    heart = [(q, sc["half_radial"] - sc["hold"], sc["half_size"] - sc["hold"])
             for tag in HEART_TAGS for r in _rows(f"{tag}/*/heart-scope.json") for q, sc in r["scores"].items()]
    embryo = [r["scores"]["fgw_shape0"] - r["scores"]["fgw"] for r in _rows("*.json", EMBRYO)]
    return {"t1": t1, "heart": heart, "embryo": embryo}


def weights(truths: list[str], own: str, w_own: float | None) -> np.ndarray:
    """Per-row weights summing to 1: each truth gets its share (w_own on `own`, the rest even; None: all even), split
    evenly over that truth's rows."""
    names = sorted(set(truths))
    share = ({q: 1 / len(names) for q in names} if w_own is None else
             {q: (w_own if q == own else (1 - w_own) / (len(names) - 1)) for q in names})
    counts = {q: truths.count(q) for q in names}
    return np.array([share[q] / counts[q] for q in truths])


def gains(d: dict, w_own: float | None) -> dict:
    t1q, t1d = zip(*d["t1"], strict=True)
    hq, htrim, hsize = zip(*d["heart"], strict=True)
    wt, wh = weights(list(t1q), "bown_cown", w_own), weights(list(hq), "own", w_own)
    we = np.full(len(d["embryo"]), 1 / len(d["embryo"]))
    a = np.asarray(t1d)
    sets = {"A": (np.asarray(htrim), np.asarray(d["embryo"])), "B": (np.asarray(hsize), np.zeros(1)),
            "C": (np.zeros(1), np.zeros(1))}
    out = {}
    for name, (b, c) in sets.items():
        wb = wh if len(b) == len(wh) else np.ones(1)
        wc = we if len(c) == len(we) else np.ones(1)
        board = (wt @ np.maximum(0, a) + 0.5 * (wb @ np.maximum(0, b) + wc @ np.maximum(0, c)))
        total = a[:, None, None] + 0.5 * (b[None, :, None] + c[None, None, :])
        per_set = np.einsum("i,j,k,ijk->", wt, wb, wc, np.maximum(0, total))
        out[name] = {"per_board": float(board), "per_set": float(per_set)}
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    d = deltas()
    res = {"rows": {"t1": len(d["t1"]), "heart": len(d["heart"]), "embryo": len(d["embryo"])},
           "weights": {}}
    lines = ["Expected gain of the s2 set over s1 alone, on the total score (T1 + mean of the T2 boards; T3 left out)",
             f"{'weight on the released-like truth':>34s}  " + "  ".join(f"{s}: per board / per set" for s in "ABC")]
    for label, w in (("even", None), ("1/3", 1 / 3), ("0.6", 0.6), ("1", 1.0)):
        g = gains(d, w)
        res["weights"][label] = g
        lines.append(f"{label:>34s}  " + "  ".join(f"{g[s]['per_board']:+6.2f} / {g[s]['per_set']:+6.2f}    "
                                                    for s in "ABC"))
    print("\n".join(lines))
    if a.out:
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "bundle_sim.json").write_text(json.dumps(res, indent=1))
        (out / "bundle_sim.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
