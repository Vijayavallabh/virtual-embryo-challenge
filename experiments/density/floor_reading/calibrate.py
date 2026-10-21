# ruff: noqa: E501  (long report lines kept as the decision panel ran this script on 2026-10-21)
"""Noise calibration for the floor-reading Monte Carlo (read-only on the repo).

Derives per-seed SDs of paired score differences between count levels from:
  1. experiments/density/review/results/analysis.txt and analysis_big.txt: 95% CI half-widths of paired means
     (t quantile, n-1 df: review_analyse.py ci()), so SD = h * sqrt(n) / t_{0.975, n-1}.
  2. experiments/density/results/results_long.csv (and v2): 5-seed per-seed rows for random and refstage at
     0.25x-8x truth_cells, from which SDs of every pairwise paired difference are computed directly.
Only pseudo-target (released training stage) numbers are used. Nothing about validation or test boards is read.

Run from the repo root:  uv run --no-sync python -I <this file>
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path.cwd()
OUT = Path(__file__).resolve().parent
DENS = REPO / "experiments" / "density"


def sd_from_halfwidth(h: float, n: int) -> float:
    return h * math.sqrt(n) / stats.t.ppf(0.975, n - 1)


def parse_analysis(path: Path):
    """Rows like 'shuffle50_2x - shuffle50_1x: heart_E9.5: +0.59 ±0.19 [0/20] | ... || 3 orig targets pooled: ...'."""
    out = []
    pat_cell = re.compile(r"([A-Za-z0-9_.]+): ([+-]\d+\.\d+) ±(\d+\.\d+) \[(\d+)/(\d+)\]")
    for line in path.read_text().splitlines():
        if " - " not in line or "±" not in line:
            continue
        head, _, rest = line.partition(": ")
        a, _, b = head.partition(" - ")
        cells = rest.split(" || ")
        for part in cells[0].split(" | "):
            m = pat_cell.search(part)
            if m:
                tgt, mean, h, k, n = m.group(1), float(m.group(2)), float(m.group(3)), int(m.group(4)), int(m.group(5))
                out.append(dict(a=a.strip(), b=b.strip(), target=tgt, mean=mean, h=h, n=n, neg=k,
                                sd=sd_from_halfwidth(h, n)))
        if len(cells) > 1:
            m = re.search(r"pooled: ([+-]\d+\.\d+) ±(\d+\.\d+) \[(\d+)/(\d+)\]", cells[1])
            if m:
                mean, h, k, n = float(m.group(1)), float(m.group(2)), int(m.group(3)), int(m.group(4))
                out.append(dict(a=a.strip(), b=b.strip(), target="pooled3", mean=mean, h=h, n=n, neg=k,
                                sd=sd_from_halfwidth(h, n)))
    return out


def ratio_of(cond: str) -> float | None:
    m = re.match(r"^[a-z0-9]+_(\d+(?:\.\d+)?)x$", cond)
    return float(m.group(1)) if m else None


def csv_pairs(path: Path, patterns=("random", "refstage", "shuffle50")):
    rows = list(csv.DictReader(path.open()))
    data = {}
    for r in rows:
        pat = r["condition"].split("_")[0]
        x = ratio_of(r["condition"])
        if pat not in patterns or x is None:
            continue
        data.setdefault((r["target"], pat), {}).setdefault(x, {})[int(r["seed"])] = float(r["score"])
    res = []
    for (tgt, pat), lv in sorted(data.items()):
        xs = sorted(lv)
        for x1, x2 in itertools.combinations(xs, 2):
            seeds = sorted(set(lv[x1]) & set(lv[x2]))
            d = np.array([lv[x2][s] - lv[x1][s] for s in seeds])
            if len(d) >= 3:
                res.append(dict(target=tgt, pattern=pat, x1=x1, x2=x2, sep=math.log2(x2 / x1), n=len(d),
                                mean=float(d.mean()), sd=float(d.std(ddof=1)), both_ge1=(x1 >= 1)))
        # level SDs after removing seed effects (two-way residuals), for the correlation structure
    return res, data


def main():
    report = {}
    lines = []
    w = lines.append

    # 1. review CI half-widths -> per-seed SDs of paired differences (shuffle50 is the rule's condition)
    rev = parse_analysis(DENS / "review" / "results" / "analysis.txt")
    big = parse_analysis(DENS / "review" / "results" / "analysis_big.txt")
    w("== per-seed SD of paired differences from review CI half-widths (SD = h*sqrt(n)/t_{.975,n-1})")
    for name, rr in (("analysis.txt (20 seeds)", rev), ("analysis_big.txt (10 seeds)", big)):
        w(f"-- {name}")
        for r in rr:
            if r["target"] == "pooled3":
                continue
            w(f"  {r['a']:>20s} - {r['b']:<14s} {r['target']:<12s} mean {r['mean']:+.2f}  h {r['h']:.2f}  n {r['n']:2d}"
              f"  -> SD {r['sd']:.3f}   [neg {r['neg']}/{r['n']}]")
    sh = [r for r in rev + big if r["a"].startswith("shuffle50") and r["b"].startswith("shuffle50")
          and r["target"] != "pooled3"]
    sh_sd = np.array([r["sd"] for r in sh])
    w(f"shuffle50 count-pair SDs (all targets, both files): n={len(sh_sd)} min {sh_sd.min():.3f} "
      f"q25 {np.quantile(sh_sd, .25):.3f} median {np.median(sh_sd):.3f} q75 {np.quantile(sh_sd, .75):.3f} "
      f"max {sh_sd.max():.3f}")
    by_t = {}
    for r in sh:
        by_t.setdefault(r["target"] + (" (big)" if r in big else ""), []).append(r["sd"])
    for t, v in by_t.items():
        w(f"  {t:<20s} median SD {np.median(v):.3f}  range {min(v):.3f}-{max(v):.3f}  (k={len(v)})")
    report["shuffle50_pair_sd"] = dict(n=len(sh_sd), min=float(sh_sd.min()), median=float(np.median(sh_sd)),
                                       q25=float(np.quantile(sh_sd, .25)), q75=float(np.quantile(sh_sd, .75)),
                                       max=float(sh_sd.max()),
                                       by_target={k: float(np.median(v)) for k, v in by_t.items()})
    # separation dependence for shuffle50 (20-seed file): 1 doubling vs 2 vs 3
    sep = {}
    for r in rev:
        if not (r["a"].startswith("shuffle50") and r["target"] != "pooled3"):
            continue
        xa, xb = ratio_of(r["a"]), ratio_of(r["b"])
        if xa and xb:
            sep.setdefault(round(abs(math.log2(xa / xb)), 2), []).append(r["sd"])
    w("shuffle50 (20 seeds) SD by separation in doublings: " +
      ", ".join(f"{k:g}: median {np.median(v):.3f} (k={len(v)})" for k, v in sorted(sep.items())))
    report["shuffle50_sd_by_sep"] = {str(k): float(np.median(v)) for k, v in sorted(sep.items())}

    # 2. 5-seed per-seed CSVs: SD of every pairwise paired difference (random, refstage), by separation
    for tag, path in (("v1", DENS / "results" / "results_long.csv"), ("v2", DENS / "results" / "v2" / "results_long.csv")):
        res, data = csv_pairs(path)
        w(f"\n== {tag}: {path.relative_to(REPO)}: per-seed paired-difference SDs, 5 seeds")
        for pat in ("random", "refstage"):
            for above in (True, False):
                sub = [r for r in res if r["pattern"] == pat and r["both_ge1"] == above]
                if not sub:
                    continue
                by = {}
                for r in sub:
                    by.setdefault(round(r["sep"], 2), []).append(r["sd"])
                w(f"  {pat:<9s} {'both levels >=1x' if above else 'lower level <1x'}: " +
                  ", ".join(f"sep {k:g}: median SD {np.median(v):.3f} (k={len(v)})" for k, v in sorted(by.items())))
                report[f"{tag}_{pat}_{'ge1' if above else 'lt1'}_sd_by_sep"] = {
                    str(k): float(np.median(v)) for k, v in sorted(by.items())}
        if tag == "v1":
            # level-noise correlation: two-way residuals per (target, pattern) over levels >= 1x
            cors = []
            for (tgt, pat), lv in data.items():
                xs = [x for x in sorted(lv) if x >= 1]
                seeds = sorted(set.intersection(*(set(lv[x]) for x in xs)))
                M = np.array([[lv[x][s] for x in xs] for s in seeds])  # seeds x levels
                R = M - M.mean(0, keepdims=True) - M.mean(1, keepdims=True) + M.mean()
                C = np.corrcoef(R.T)
                for i, j in itertools.combinations(range(len(xs)), 2):
                    cors.append((tgt, pat, xs[i], xs[j], C[i, j]))
            w("  two-way residual correlations between levels >=1x are biased toward -1/(L-1) with 5 seeds; "
              "not used. Correlation structure is instead read from SD growth with separation above.")
        # compare v1/v2 identical rows for random/refstage
    v1 = {(r["target"], r["seed"], r["condition"]): r["score"] for r in csv.DictReader((DENS / "results" / "results_long.csv").open())}
    v2 = {(r["target"], r["seed"], r["condition"]): r["score"] for r in csv.DictReader((DENS / "results" / "v2" / "results_long.csv").open())}
    common = [k for k in v1 if k in v2 and k[2].split("_")[0] in ("random", "refstage")]
    same = sum(abs(float(v1[k]) - float(v2[k])) < 1e-9 for k in common)
    w(f"\nv1 vs v2 random/refstage rows identical: {same}/{len(common)}")

    # 3. typical true gains offline (shuffle50 means) for the 'more is better' family calibration
    w("\n== shuffle50 paired mean gains offline (for the increasing family's slope range)")
    for r in rev + big:
        if r["a"].startswith("shuffle50") and r["b"] == "shuffle50_1x":
            xa = ratio_of(r["a"])
            w(f"  {r['a']:<14s} vs 1x {r['target']:<12s} {r['mean']:+.2f}  -> per doubling {r['mean'] / math.log2(xa):+.2f}")

    (OUT / "calibration.txt").write_text("\n".join(lines) + "\n")
    (OUT / "calibration.json").write_text(json.dumps(report, indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
