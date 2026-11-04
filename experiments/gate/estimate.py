#!/usr/bin/env python3
"""Our estimated leaderboard total for the captain's upload gate, from committed results only.

Aggregation only, not a try. Every input is a committed practice-task result (JSON or CSV under experiments/*/results/),
read through one door (`Door`) that refuses boards/, panels/, docs/research/, data/, runs/ and any index.json; the one
other file is results/upload_gate.json, for the captain's threshold. Nothing here reads a target, a test-board file or
a board endpoint, runs a candidate, a builder or a scorer, or computes a new score. From the files it reads it uses
only labels, seeds and board scores; no cell count, anchor or target field is used or printed. Every number is a mean,
a maximum or a percentile of board scores already committed. No leaderboard figure enters it, as a constant or
otherwise (docs/PLAN.md section 8; docs/DATA_POLICY.md 2026-11-05).

The two readings of "our estimated score in leaderboard" (docs/PLAN.md:605-636):
  A  every board at the better of its two uploads, since a board shows a team's best;
  B  the first round alone: T1 s1, heart s1, embryo s1 and T3's s2 (on T3 the first upload is s2).
Overall = T1 + (heart + embryo) / 2 + T3.

The procedure follows the evidence decision log 2026-11-04 (c) names (docs/PLAN.md:989-993), with the simplest rule
where (c) states no step (experiments/gate/README.md); it does not reproduce (c)'s figures:
  T1      the oracle-base grid of experiments/s2_hedge/ under its uniform prior (s1 bown_cown, s2 bown_c0.6), plus
          M1's paired gain over B0 on the same backtest (diag_v1), whose B0 is the harness's (d2_v1) seed for seed;
  heart   the heart-scope payoffs under their uniform prior (s1 hold, s2 half_radial) on the two transitions the
          harness counts, each moved so that s1 on the released truth scores what the entry scored on that backtest
          (d2_v1's b0, without the unscored size): the payoffs' oracle core size gets no credit;
  embryo  the two counted backtests (emb_shape_v1; s1 fgw, s2 fgw_shape0);
  T3      s1's two stand-ins on the Gata4 replicates, evenly, against s2's stand-in (gata4_records.csv); with K2a
          (decision log 2026-11-04 (d)) the cross-knockout half of s1's weight is shared evenly by its two directions,
          stand-in B (the Mab21l2 factor scored on Gata4) and K2a (s1's Gata4 factor scored on Mab21l2, against our
          WT copy on that fold for s2).
A board's outcomes are committed per-seed scores, one per scenario (truth, transition, backtest, stand-in, replicate)
and seed, weighted by the procedure. A reading is their weighted mean. Its band is the 10th and 90th percentiles:
exact on each board's weighted outcomes, and for Overall a fixed-seed Monte Carlo that draws one outcome per board,
the boards independent.

Usage, from the repository root:
  uv run --frozen --no-sync python -I -B experiments/gate/estimate.py [--json PATH]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]

# ------------------------------------------------------------------------------ the committed outlook, compared only
# Decision log 2026-11-04 (c), docs/PLAN.md:989-991: "about 184 overall, likely 170-198 (T1 about 63, T3 about 65, heart
# about 52, embryo about 60, each board at the better of its two uploads)"; results/upload_gate.json:5-9 and :11. No
# step below is fitted to these figures; they are printed beside the readings and differenced.
OUTLOOK = {"T1": 63, "heart": 52, "embryo": 60, "T3": 65}  # docs/PLAN.md:990
OUTLOOK_TOTAL = 184  # docs/PLAN.md:990
OUTLOOK_RANGE = (170, 198)  # docs/PLAN.md:990
AGREE_WITHIN = 1  # points per board for "reproduces the outlook" (chosen: experiments/gate/README.md, The procedure)

# ------------------------------------------------------------------------------------------- readings and Overall
BOARDS = ("T1", "heart", "embryo", "T3")
OVERALL_WEIGHT = {"T1": 1.0, "heart": 0.5, "embryo": 0.5, "T3": 1.0}  # docs/PLAN.md:609; docs/CHALLENGE.md:109-112
FIRST_ROUND = {"T1": "s1", "heart": "s1", "embryo": "s1", "T3": "s2"}  # docs/PLAN.md:607, :636; docs/S2_RULES.md:282
PRIOR = "uniform"  # every payoff's first stated prior: experiments/s2_hedge/payoff.py:483; docs/S2_RULES.md:95, :144
BAND = (10, 90)  # the band's percentiles (chosen: experiments/gate/README.md, The procedure)
N_DRAWS = 200_000  # Monte Carlo draws behind the bands
RNG_SEED = 20261105  # fixed so the bands reproduce (the date the estimator was written)
TOL = 1e-9  # largest gap allowed when a committed summary is recomputed from its per-seed files

# ------------------------------------------------------------------------------------------------------------- T1
# s1 is M1 at E10.5's own chamber mix, s2 the same at cardiomyocytes 0.60 of heart cells (docs/S2_RULES.md:85). That bet
# (docs/PLAN.md:991-992: "E10.5's chamber mix with the s2 hedge") is priced by the oracle-base grid
# (experiments/s2_hedge/payoff.py:28-32): B0 built at the truth's own chamber mix with the two shares set.
T1_TASK = "T1/E9.5->E10.5:chamber"
T1_GRID_RUNS = "experiments/s2_hedge/results/per_seed/t1-grid"  # experiments/s2_hedge/README.md:220-224
T1_GRID_PAYOFF = "experiments/s2_hedge/results/t1_grid/payoff.json"  # mean[truth][pred]; checked against the runs
T1_GRID_HEDGE = "experiments/s2_hedge/results/t1_grid/hedge_bown_cown.json"  # its priors; checked against the runs
T1_S1 = "bown_cown"  # experiments/s2_hedge/README.md:54 (s1's analogue)
T1_S2 = "bown_c0.6"  # docs/S2_RULES.md:85 (cardiomyocytes 0.60 of heart cells, blood at its own share)
# The backtest-to-test step: the grid is B0 and the entries are M1 (B0 plus Feng's direction), so M1's paired gain over
# B0 on the chamber-cut backtest is added to every T1 outcome (experiments/harness/TRIES.md:99;
# experiments/t1_m1/results/diag_v1/README.md); diag_v1's B0 must be the harness's B0 (d2_v1) seed for seed.
T1_DIAG = "experiments/t1_m1/results/diag_v1/T1_E9.5_E10.5_chamber.json"
T1_D2 = "experiments/harness/results/d2_v1/tasks/T1_E9.5_E10.5_chamber.json"
T1_B0, T1_M1 = "b0", "m1_feng"

# ---------------------------------------------------------------------------------------------------------- heart
# s1 is H1, s2 the half trim (docs/S2_RULES.md:118-121): hold and half_radial in the payoffs
# (experiments/s2_hedge/payoff.py:21-27; experiments/s2_hedge/README.md:120-122). The payoffs scale every prediction to
# the truth's core size (payoff.py:25-27, an oracle size); the entry's pre-registered size is scored by no backtest
# (experiments/harness/TRIES.md:133-135, A8). So each transition is moved by the entry's score on its counted backtest
# (d2_v1's b0, TRIES.md:115-116) minus hold's on the released truth, and the payoffs add only what the truths and the s2
# change. E8.25_late -> E8.5 has no counted backtest and stays out (experiments/s2_hedge/README.md:166-168).
HEART = (  # (per-seed payoff runs, their hedge file, the counted backtest of the same transition)
    ("experiments/s2_hedge/results/per_seed/heart-scope", "experiments/s2_hedge/results/heart_scope/hedge_hold.json",
     "experiments/harness/results/d2_v1/tasks/T2h_E9.5_E10.5.json"),
    ("experiments/s2_hedge/results/per_seed/heart-scope-e875",
     "experiments/s2_hedge/results/heart_scope_e875/hedge_hold.json",
     "experiments/harness/results/d2_v1/tasks/T2h_E8.75_E10.5.json"),
)
HEART_S1, HEART_S2 = "hold", "half_radial"
HEART_RELEASED = "own"  # experiments/s2_hedge/payoff.py:21 (E10.5 as released)
HEART_ENTRY = "b0"  # experiments/harness/TRIES.md:115 (H1*: b0, H1 without the pre-registered size)

# --------------------------------------------------------------------------------------------------------- embryo
# s1 fgw, s2 fgw_shape0 (docs/S2_RULES.md:239-255) on the two counted backtests (experiments/harness/TRIES.md:117-120);
# emb_shape_v1's fgw must be emb_fgw_v1's seed for seed (experiments/s2_hedge/README.md:208-209).
EMBRYO = (  # (the s2 rule's run, the s1 run, on the same backtest)
    ("experiments/harness/results/emb_shape_v1/tasks/T2e_E7.25_E8.0_E7.5.json",
     "experiments/harness/results/emb_fgw_v1/tasks/T2e_E7.25_E8.0_E7.5.json"),
    ("experiments/harness/results/emb_shape_v1/tasks/T2e_E6.75_E7.5_E7.25.json",
     "experiments/harness/results/emb_fgw_v1/tasks/T2e_E6.75_E7.5_E7.25.json"),
)
EMBRYO_S1, EMBRYO_S2 = "fgw", "fgw_shape0"

# ------------------------------------------------------------------------------------------------------------- T3
# No backtest runs the T3 entries; their stand-ins are rows of the pre-registered decomposition on the two Gata4
# replicates (experiments/harness/TRIES.md:137-146). s1's bet, the knockout-depth term (docs/PLAN.md:992), is in
# stand-in A and not in B; the two count evenly, a chosen rule, since (c) states no weights (experiments/gate/README.md,
# The procedure).
T3_RECORDS = "experiments/t3_decomp/results/gata4/gata4_records.csv"  # per seed; only target, seed, variant, t3 read
T3_SUMMARY = "experiments/t3_decomp/results/gata4/gata4_summary.csv"  # seed means; checked against the records
T3_REPLICATES = ("rep1", "rep2")
T3_STANDINS = {"A": "xrep:F+depth", "B": "transfer:mab21l2_all"}  # experiments/harness/TRIES.md:145, :144
T3_S2 = "blind:G2_bcat_central_d50"  # experiments/harness/TRIES.md:146
T3_K2 = "experiments/harness/results/t3_k2_v1/tasks/T3_Gata4_Mab21l2.json"  # docs/PLAN.md:1022-1026
T3_K2_TASK = "T3/Gata4->Mab21l2"
T3_K2A = "t3_factor_xko"  # experiments/harness/TRIES.md:258 (row K2a)
T3_K2_S2 = "copy"  # our WT copy on K2a's fold, standing in there for T3's near-WT s2 (chosen; README, The procedure)
# K2a, included because docs/PLAN.md:1091-1092 asks for a symmetric restatement; the split is chosen: the cross-knockout
# half of s1's weight is split evenly between its two directions, B (Mab21l2 -> Gata4) and K2a (Gata4 -> Mab21l2), so a
# K2a above B would have raised T3 exactly as this one lowers it. K2a's fold has no s2 stand-in, so our WT copy on that
# fold stands in for s2, paired by seed: both uploads then meet the same truth, as a board's better of two does.
T3_WEIGHTS_AS_OF_C = {"A": 0.5, "B": 0.5}
T3_WEIGHTS_K2A = {"A": 0.5, "B": 0.25, "K2a": 0.25}

GATE = "results/upload_gate.json"  # the captain's threshold, read and never written (AGENTS.md:95-104)

# ---------------------------------------------------------------------------------------------------- the one door
READABLE = ("experiments/s2_hedge/results/", "experiments/harness/results/", "experiments/t1_m1/results/diag_v1/",
            "experiments/t3_decomp/results/gata4/", GATE)
REFUSED = ("boards/", "panels/", "docs/research/", "data/", "runs/")


class Door:
    """The one way in: committed results under READABLE, never a board file; every file read is hashed."""

    def __init__(self, root: Path = REPO):
        self.root = Path(root).resolve()
        self.read: dict[str, str] = {}

    def _allowed(self, rel: str) -> Path:
        readable = any(rel == p or (p.endswith("/") and rel.startswith(p)) for p in READABLE)
        if (rel.startswith(REFUSED) or "index.json" in rel or ".." in Path(rel).parts or Path(rel).is_absolute()
                or not readable):
            raise SystemExit(f"estimate: refusing {rel}: not a committed result this estimator reads")
        path = (self.root / rel).resolve()
        if self.root not in path.parents:
            raise SystemExit(f"estimate: refusing {rel}: outside the repository")
        real = path.relative_to(self.root).as_posix()
        if real != rel.rstrip("/") and (real.startswith(REFUSED) or "index.json" in real or not any(
                real == p.rstrip("/") or (p.endswith("/") and real.startswith(p)) for p in READABLE)):
            raise SystemExit(f"estimate: refusing {rel}: it resolves to {real}")
        return path

    def path(self, rel: str) -> Path:
        path = self._allowed(rel)
        if not path.is_file():
            raise SystemExit(f"estimate: missing committed input {rel}")
        self.read[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    def json(self, rel: str):
        return json.loads(self.path(rel).read_text())

    def runs(self, rel_dir: str) -> list[str]:
        """A payoff's per-seed run files (one folder per part), sorted."""
        files = sorted(self._allowed(rel_dir.rstrip("/") + "/").glob("*/*.json"))
        if not files:
            raise SystemExit(f"estimate: missing committed input {rel_dir}/*/*.json")
        return [f.relative_to(self.root).as_posix() for f in files]


def required() -> list[str]:
    """Every committed input: files, and the per-seed run folders (ending in /)."""
    files = [T1_GRID_PAYOFF, T1_GRID_HEDGE, T1_DIAG, T1_D2, T3_RECORDS, T3_SUMMARY, T3_K2, GATE]
    files += [x for h in HEART for x in h[1:]] + [x for e in EMBRYO for x in e]
    return files + [T1_GRID_RUNS + "/"] + [h[0] + "/" for h in HEART]


def missing_inputs(root: Path = REPO) -> list[str]:
    root = Path(root)
    out = []
    for rel in required():
        if rel.endswith("/"):
            if not any((root / rel).glob("*/*.json")):
                out.append(rel + "*/*.json")
        elif not (root / rel).is_file():
            out.append(rel)
    return out


# -------------------------------------------------------------------------------------------------------- loaders
def load_runs(door: Door, rel_dir: str) -> dict:
    """A payoff's per-seed runs (experiments/s2_hedge/payoff.py), merged over parts as its hedge() merges them:
    {"task", "truths", "preds", "seeds", "scores": {truth: {pred or "floor": per-seed array}}}. Only the labels and
    each row's seed and scores are read."""
    merged: dict[int, dict] = {}
    head = None
    for rel in door.runs(rel_dir):
        d = door.json(rel)
        try:
            this = (d["task"], tuple(d["truth_values"]), tuple(d["pred_values"]))
            rows = [(int(r["seed"]), r["scores"]) for r in d["rows"]]
        except (KeyError, TypeError, ValueError) as e:
            raise SystemExit(f"estimate: {rel}: missing or malformed field {e}") from None
        if head is None:
            head = this
        elif this != head:
            raise SystemExit(f"estimate: {rel} is not a part of the same payoff as the rest of {rel_dir}")
        for seed, scores in rows:
            merged.setdefault(seed, {}).update(scores)
    task, truths, preds = head
    seeds = sorted(merged)
    for k in seeds:
        for t in truths:
            if t not in merged[k] or any(p not in merged[k][t] for p in (*preds, "floor")):
                raise SystemExit(f"estimate: {rel_dir}: seed {k} lacks a score for truth {t}")
    try:
        scores = {t: {p: np.array([merged[k][t][p] for k in seeds], dtype=float) for p in (*preds, "floor")}
                  for t in truths}
    except (KeyError, TypeError, ValueError) as e:
        raise SystemExit(f"estimate: {rel_dir}: missing or malformed score {e}") from None
    return {"task": task, "truths": list(truths), "preds": list(preds), "seeds": seeds, "scores": scores}


def load_task(door: Door, rel: str, cands: tuple[str, ...]) -> dict:
    """A committed harness or diagnostic task file: the task's key and each named candidate's board score per seed
    (seeds sorted). Only the task key and each row's seed and scores are used."""
    d = door.json(rel)
    try:
        rows = sorted(d["rows"], key=lambda r: r["seed"])
        key = d["task"]["key"] if isinstance(d["task"], dict) else d["task"]
        out = {"key": key, "seeds": [int(r["seed"]) for r in rows]}
        for c in cands:
            out[c] = np.array([float(r["scores"][c]) for r in rows])
    except (KeyError, TypeError, ValueError) as e:
        raise SystemExit(f"estimate: {rel}: missing or malformed field {e}") from None
    return out


def load_t3(door: Door, variants: tuple[str, ...], checks: list[str]) -> dict:
    """Per-seed board scores (column t3) of the named decomposition variants on each Gata4 replicate, checked against
    the summary's seed means. Only the columns target, seed, variant and t3 are used."""
    per: dict = {}
    means: dict = {}
    for rel in (T3_RECORDS, T3_SUMMARY):
        try:
            with open(door.path(rel), newline="") as fh:
                for r in csv.DictReader(fh):
                    if r["variant"] not in variants:
                        continue
                    if rel == T3_RECORDS:
                        per.setdefault(r["target"], {}).setdefault(r["variant"], {})[int(r["seed"])] = float(r["t3"])
                    else:
                        means[r["target"], r["variant"]] = float(r["t3"])
        except (KeyError, TypeError, ValueError) as e:
            raise SystemExit(f"estimate: {rel}: missing or malformed field {e}") from None
    out = {}
    for target in T3_REPLICATES:
        by_var = per.get(target, {})
        lacking = [v for v in variants if v not in by_var or (target, v) not in means]
        if lacking:
            raise SystemExit(f"estimate: {T3_RECORDS} or {T3_SUMMARY}: {target} lacks {lacking}")
        seeds = sorted(by_var[variants[0]])
        if any(sorted(by_var[v]) != seeds for v in variants):
            raise SystemExit(f"estimate: {T3_RECORDS}: {target}'s variants have different seeds")
        out[target] = {v: np.array([by_var[v][s] for s in seeds]) for v in variants}
        same([out[target][v].mean() for v in variants], [means[target, v] for v in variants],
             f"T3 {target}: the records' seed means are gata4_summary.csv's t3 ({len(variants)} variants)", checks)
    return out


# --------------------------------------------------------------------------------------------------- the procedure
def same(x, y, what: str, checks: list[str], tol: float = TOL) -> None:
    """A self-check: x and y agree within tol everywhere, or the estimate stops."""
    gap = float(np.max(np.abs(np.asarray(x, dtype=float) - np.asarray(y, dtype=float))))
    if not gap <= tol:
        raise SystemExit(f"estimate: check failed: {what} (largest gap {gap:.3g})")
    checks.append(f"{what} (largest gap {gap:.1e})")


def pack(a: list, b: list, w: list) -> dict:
    w_ = np.concatenate(w)
    return {"a": np.concatenate(a), "b": np.concatenate(b), "w": w_ / w_.sum()}


def grid_outcomes(s1: dict, s2: dict, prior: dict[str, float], shift: float = 0.0) -> dict:
    """One outcome per truth and seed of a payoff, weighted prior(truth) / seeds: "a" is the better of s1 and s2 in that
    seed (the board shows the better upload), "b" is s1 alone, both moved by `shift`. s1 and s2 map each truth to its
    per-seed scores, paired by seed."""
    a, b, w = [], [], []
    for t, p in prior.items():
        x1, x2 = np.asarray(s1[t], dtype=float), np.asarray(s2[t], dtype=float)
        if x1.shape != x2.shape:
            raise ValueError(f"truth {t}: s1 and s2 are not paired by seed")
        a.append(np.maximum(x1, x2) + shift)
        b.append(x1 + shift)
        w.append(np.full(x1.size, p / x1.size))
    return pack(a, b, w)


def mixture(parts: list[tuple[float, dict]]) -> dict:
    """Outcome sets mixed by weight (the weights are normalised)."""
    tot = sum(p for p, _ in parts)
    return pack([o["a"] for _, o in parts], [o["b"] for _, o in parts], [o["w"] * (p / tot) for p, o in parts])


def standin_outcomes(s1, s2, paired: bool) -> dict:
    """T3 outcomes of one s1 stand-in against s2's stand-in: "a" is the better of the two uploads, "b" is s2 alone (T3's
    first upload). Paired by seed when both were scored on the same replicate by the same scorer; otherwise every seed
    of one meets every seed of the other."""
    s1, s2 = np.asarray(s1, dtype=float), np.asarray(s2, dtype=float)
    if paired:
        if s1.shape != s2.shape:
            raise ValueError("paired stand-ins need the same seeds")
        x1, x2 = s1, s2
    else:
        x1, x2 = np.repeat(s1, s2.size), np.tile(s2, s1.size)
    return {"a": np.maximum(x1, x2), "b": x2, "w": np.full(x1.size, 1.0 / x1.size)}


def t3_outcomes(standins: dict, s2: dict, weights: dict[str, float]) -> dict:
    """standins: name -> {replicate: per-seed scores} for a stand-in on the Gata4 replicates (paired with s2's stand-in
    on the same replicate), or a pair (per-seed scores, the s2 stand-in's per-seed scores on the same fold) for one
    scored on another fold (K2a, against our WT copy there); s2: {replicate: per-seed scores}. The replicates count
    evenly."""
    parts = []
    for name, wt in weights.items():
        x = standins[name]
        if isinstance(x, dict):
            reps = [standin_outcomes(x[r], s2[r], paired=True) for r in s2]
        elif isinstance(x, tuple) and len(x) == 2:
            reps = [standin_outcomes(x[0], x[1], paired=True)]
        else:
            raise TypeError(f"stand-in {name}: a dict by replicate or a (scores, s2 scores) pair on its own fold")
        parts.append((wt, mixture([(1.0, o) for o in reps])))
    return mixture(parts)


def reading(o: dict, key: str) -> float:
    return float(np.dot(o[key], o["w"]))


def overall(values: dict[str, float]) -> float:
    return sum(OVERALL_WEIGHT[b] * values[b] for b in BOARDS)


def readings(boards: dict[str, dict]) -> dict:
    out = {}
    for name, key in (("A", "a"), ("B", "b")):
        r = {b: reading(boards[b], key) for b in BOARDS}
        out[name] = {**r, "Overall": overall(r)}
    return out


def weighted_percentile(x, w, q: float) -> float:
    """The q-th percentile of outcomes x with weights w, exactly: the smallest outcome whose cumulative weight reaches
    q; where the cumulative weight meets q exactly on an outcome, the midpoint of that outcome and the next."""
    order = np.argsort(x, kind="stable")
    xs, cw = np.asarray(x, dtype=float)[order], np.cumsum(np.asarray(w, dtype=float)[order])
    cw /= cw[-1]
    target = q / 100
    i = int(np.searchsorted(cw, target - TOL))
    if abs(cw[i] - target) <= TOL and i + 1 < xs.size:
        return float((xs[i] + xs[i + 1]) / 2)
    return float(xs[i])


def bands(boards: dict[str, dict], n: int = N_DRAWS, seed: int = RNG_SEED) -> dict:
    """The BAND percentiles: exact on each board's weighted outcomes, and of Overall from n draws of one outcome per
    board (the boards independent; A and B read the same drawn scenario)."""
    rng = np.random.default_rng(seed)
    out: dict = {"A": {}, "B": {}}
    tot = {"A": np.zeros(n), "B": np.zeros(n)}
    for b in BOARDS:
        o = boards[b]
        idx = rng.choice(o["w"].size, size=n, p=o["w"])
        for name, key in (("A", "a"), ("B", "b")):
            out[name][b] = [weighted_percentile(o[key], o["w"], q) for q in BAND]
            tot[name] += OVERALL_WEIGHT[b] * o[key][idx]
    for name in out:
        out[name]["Overall"] = [float(v) for v in np.percentile(tot[name], BAND)]
    return out


def hedge_value(grid: dict, s1: str, s2: str, prior: dict[str, float]) -> float:
    """experiments/s2_hedge/payoff.py hedge(): the seed mean of the prior-weighted better of s1 and s2."""
    sc = grid["scores"]
    return reading(grid_outcomes({t: sc[t][s1] for t in prior}, {t: sc[t][s2] for t in prior}, prior), "a")


def check_hedge(grid: dict, hedge: dict, s1: str, s2: str, what: str, checks: list[str]) -> None:
    try:
        if hedge["s1"] != s1 or [int(k) for k in hedge["seeds"]] != grid["seeds"]:
            raise SystemExit(f"estimate: {what}: the hedge file is not s1 = {s1} on the seeds of the per-seed runs")
        wants = {cand: [hedge["candidates"][cand][name]["mean"] for name in hedge["priors"]] for cand in (s1, s2)}
    except (KeyError, TypeError, ValueError) as e:
        raise SystemExit(f"estimate: {what}: missing or malformed field {e}") from None
    for cand in (s1, s2):
        got = [hedge_value(grid, s1, cand, prior) for prior in hedge["priors"].values()]
        want = wants[cand]
        value = f"{s1} alone" if cand == s1 else f"the better of {s1} and {cand}"
        same(got, want, f"{what}: the per-seed runs give the hedge file's {value} under {len(want)} stated priors",
             checks)


# ----------------------------------------------------------------------------------------------------------- boards
def board_t1(door: Door, checks: list[str]) -> tuple[dict, dict, dict]:
    grid = load_runs(door, T1_GRID_RUNS)
    payoff, hedge = door.json(T1_GRID_PAYOFF), door.json(T1_GRID_HEDGE)
    if not grid["task"] == payoff["task"] == T1_TASK or (payoff["truth_values"], payoff["pred_values"]) != (
            grid["truths"], grid["preds"]):
        raise SystemExit(f"estimate: {T1_GRID_PAYOFF} and {T1_GRID_RUNS} are not the same {T1_TASK} grid")
    truths = grid["truths"]
    for p in (T1_S1, T1_S2):
        same([grid["scores"][t][p].mean() for t in truths], [payoff["mean"][t][p] for t in truths],
             f"T1 grid: the per-seed runs give payoff.json's mean[truth][{p}] on {len(truths)} truths", checks)
    check_hedge(grid, hedge, T1_S1, T1_S2, f"T1 grid ({T1_GRID_HEDGE})", checks)
    diag, d2 = load_task(door, T1_DIAG, (T1_B0, T1_M1)), load_task(door, T1_D2, (T1_B0, "floor"))
    if not diag["key"] == d2["key"] == T1_TASK or diag["seeds"] != d2["seeds"]:
        raise SystemExit(f"estimate: {T1_DIAG} and {T1_D2} are not the same {T1_TASK} seeds")
    same(np.concatenate([grid["scores"][t]["floor"] for t in truths] + [d2["floor"]]), 50.0,
         "T1: the floor is 50 on every truth and seed of the grid and on the backtest", checks)
    same(diag[T1_B0], d2[T1_B0], "T1: diag_v1's B0 is the harness's B0 (d2_v1) seed for seed", checks)
    step = float(np.mean(diag[T1_M1] - diag[T1_B0]))
    sc = grid["scores"]

    def under(prior: dict[str, float]) -> dict:
        return grid_outcomes({t: sc[t][T1_S1] for t in prior}, {t: sc[t][T1_S2] for t in prior}, prior, shift=step)

    by_prior = {name: under(prior) for name, prior in hedge["priors"].items()}
    steps = {"m1_minus_b0": step, "grid_A": hedge_value(grid, T1_S1, T1_S2, hedge["priors"][PRIOR]),
             "grid_B": hedge_value(grid, T1_S1, T1_S1, hedge["priors"][PRIOR])}
    return by_prior[PRIOR], by_prior, steps


def board_heart(door: Door, checks: list[str]) -> tuple[dict, dict, list[dict]]:
    by_prior: dict[str, list] = {}
    priors = None
    offsets = []
    for runs_dir, hedge_rel, d2_rel in HEART:
        g, hedge = load_runs(door, runs_dir), door.json(hedge_rel)
        d2 = load_task(door, d2_rel, (HEART_ENTRY, "floor"))
        if g["task"] != d2["key"]:
            raise SystemExit(f"estimate: {runs_dir} is {g['task']}, but {d2_rel} is {d2['key']}")
        if priors is not None and hedge["priors"] != priors:
            raise SystemExit(f"estimate: {hedge_rel} states other priors than the first heart transition")
        priors = hedge["priors"]
        check_hedge(g, hedge, HEART_S1, HEART_S2, f"heart {g['task']} ({hedge_rel})", checks)
        # the floor is 50 on the released truth and the backtest (on a synthetic trim it can sit above 50 in a seed)
        same(np.concatenate([g["scores"][HEART_RELEASED]["floor"], d2["floor"]]), 50.0,
             f"heart {g['task']}: the floor is 50 on the released truth and on the backtest", checks)
        level, released = float(d2[HEART_ENTRY].mean()), float(g["scores"][HEART_RELEASED][HEART_S1].mean())
        offsets.append({"task": g["task"], "backtest_b0": level, "payoff_hold_released": released,
                        "offset": level - released})
        sc = g["scores"]
        for name, prior in priors.items():
            by_prior.setdefault(name, []).append(grid_outcomes(
                {t: sc[t][HEART_S1] for t in prior}, {t: sc[t][HEART_S2] for t in prior}, prior,
                shift=level - released))
    mixed = {name: mixture([(1.0, o) for o in parts]) for name, parts in by_prior.items()}
    return mixed[PRIOR], mixed, offsets


def board_embryo(door: Door, checks: list[str]) -> dict:
    parts = []
    for shape_rel, fgw_rel in EMBRYO:
        s = load_task(door, shape_rel, ("floor", EMBRYO_S1, EMBRYO_S2))
        f = load_task(door, fgw_rel, (EMBRYO_S1,))
        if s["key"] != f["key"] or s["seeds"] != f["seeds"]:
            raise SystemExit(f"estimate: {shape_rel} and {fgw_rel} are not the same backtest and seeds")
        same(s[EMBRYO_S1], f[EMBRYO_S1], f"embryo {s['key']}: emb_shape_v1's fgw is emb_fgw_v1's seed for seed",
             checks)
        same(s["floor"], 50.0, f"embryo {s['key']}: the floor is 50 on every seed", checks)
        parts.append((1.0, grid_outcomes({"backtest": s[EMBRYO_S1]}, {"backtest": s[EMBRYO_S2]}, {"backtest": 1.0})))
    return mixture(parts)


def board_t3(door: Door, checks: list[str]) -> tuple[dict, dict, dict]:
    rec = load_t3(door, (T3_STANDINS["A"], T3_STANDINS["B"], T3_S2), checks)
    k2 = load_task(door, T3_K2, ("floor", T3_K2A, T3_K2_S2))
    if k2["key"] != T3_K2_TASK:
        raise SystemExit(f"estimate: {T3_K2} is {k2['key']}, not {T3_K2_TASK}")
    same(k2["floor"], 50.0, f"T3 {T3_K2_TASK}: the floor is 50 on every seed", checks)
    standins = {name: {r: rec[r][v] for r in T3_REPLICATES} for name, v in T3_STANDINS.items()}
    standins["K2a"] = (k2[T3_K2A], k2[T3_K2_S2])
    s2 = {r: rec[r][T3_S2] for r in T3_REPLICATES}
    means = {name: float(np.mean([np.mean(x[r]) for r in T3_REPLICATES])) for name, x in standins.items()
             if isinstance(x, dict)}
    means.update({"K2a": float(np.mean(k2[T3_K2A])), "K2a_fold_copy": float(np.mean(k2[T3_K2_S2])),
                  "s2": float(np.mean([np.mean(s2[r]) for r in s2]))})
    return t3_outcomes(standins, s2, T3_WEIGHTS_AS_OF_C), t3_outcomes(standins, s2, T3_WEIGHTS_K2A), means


# --------------------------------------------------------------------------------------------------------- estimate
def estimate(root: Path = REPO, n_draws: int = N_DRAWS, seed: int = RNG_SEED) -> dict:
    """Both readings, per board and Overall, with their bands, the outlook comparison and every self-check."""
    missing = missing_inputs(root)
    if missing:
        raise SystemExit("estimate: missing committed inputs, so nothing is estimated:\n  " + "\n  ".join(missing))
    door = Door(root)
    checks: list[str] = []
    t1, t1_priors, t1_steps = board_t1(door, checks)
    heart, heart_priors, heart_offsets = board_heart(door, checks)
    embryo = board_embryo(door, checks)
    t3_c, t3_k, t3_means = board_t3(door, checks)
    boards_c = {"T1": t1, "heart": heart, "embryo": embryo, "T3": t3_c}
    boards_k = {**boards_c, "T3": t3_k}
    r_c, r_k = readings(boards_c), readings(boards_k)
    # each T3 scenario is one truth with an s1 and an s2 stand-in; as of (c) both scenarios are the Gata4 replicates,
    # so reading B's T3 is s2's stand-in there (with K2a, a quarter of it is our WT copy on K2a's fold)
    same(r_c["B"]["T3"], t3_means["s2"], "reading B as of (c): T3 is s2's stand-in on the Gata4 replicates", checks)
    b_c, b_k = bands(boards_c, n_draws, seed), bands(boards_k, n_draws, seed)
    diff_c = {b: r_c["A"][b] - OUTLOOK[b] for b in BOARDS}
    diff_c["Overall"] = r_c["A"]["Overall"] - OUTLOOK_TOTAL
    diff_k = {b: r_k["A"][b] - OUTLOOK[b] for b in BOARDS}
    diff_k["Overall"] = r_k["A"]["Overall"] - OUTLOOK_TOTAL
    threshold = door.json(GATE).get("threshold")
    if not isinstance(threshold, int | float):
        raise SystemExit(f"estimate: {GATE} holds no threshold")
    sens = {"T1": {n: {"A": reading(o, "a"), "B": reading(o, "b")} for n, o in t1_priors.items()},
            "heart": {n: {"A": reading(o, "a"), "B": reading(o, "b")} for n, o in heart_priors.items()}}
    return {
        "label": "our estimated leaderboard total from committed results only: aggregation, not a try; no test score "
                 "and no leaderboard figure",
        "readings": {"A": {**r_k["A"], "band": b_k["A"]}, "A_as_of_c": {**r_c["A"], "band": b_c["A"]},
                     "B": {**r_k["B"], "band": b_k["B"]}, "B_as_of_c": {**r_c["B"], "band": b_c["B"]}},
        "first_round": FIRST_ROUND,
        "outlook": {**OUTLOOK, "Overall": OUTLOOK_TOTAL, "range": list(OUTLOOK_RANGE),
                    "source": "docs/PLAN.md:989-991 (decision log 2026-11-04 (c))"},
        "A_as_of_c_minus_outlook": diff_c, "A_minus_outlook": diff_k,
        "reproduces_outlook": bool(all(abs(diff_c[b]) <= AGREE_WITHIN for b in BOARDS)),
        "agree_within": AGREE_WITHIN,
        "steps": {"prior": PRIOR, "T1": t1_steps, "heart": heart_offsets, "T3_stand_in_means": t3_means,
                  "T3_weights_as_of_c": T3_WEIGHTS_AS_OF_C, "T3_weights_with_k2a": T3_WEIGHTS_K2A},
        "sensitivity_to_stated_priors": sens,
        "gate": {"threshold": threshold, "source": f"{GATE} (the captain's upload gate)",
                 "A_above": bool(r_k["A"]["Overall"] > threshold),
                 "A_as_of_c_above": bool(r_c["A"]["Overall"] > threshold),
                 "B_above": bool(r_k["B"]["Overall"] > threshold),
                 "note": "a crossing also needs an independent recomputation, docs/PLAN.md section 8"},
        "band_percentiles": list(BAND), "draws": n_draws, "rng_seed": seed,
        "checks": checks, "inputs_sha256": dict(sorted(door.read.items())),
    }


# ----------------------------------------------------------------------------------------------------------- output
def _wrap(head: str, items: list[str], width: int = 118) -> list[str]:
    lines, cur = [], head
    for it in items:
        if len(cur) + len(it) + 2 > width and cur.strip():
            lines.append(cur.rstrip())
            cur = " " * len(head)
        cur += it + "  "
    lines.append(cur.rstrip())
    return lines


def table(res: dict) -> str:
    r = res["readings"]

    def row(label: str, d: dict, band=None, signed: bool = False) -> str:
        f = "{:+8.2f}" if signed else "{:8.2f}"
        s = f"{label:46s}" + "".join(f.format(d[b]) for b in BOARDS) + ("{:+9.2f}" if signed else "{:9.2f}").format(
            d["Overall"])
        return s + (f"  [{band[0]:.1f}, {band[1]:.1f}]" if band is not None else "")

    out = ["Our estimated leaderboard total: aggregation of committed results, not a try "
           "(experiments/gate/estimate.py)",
           f"{'':46s}" + "".join(f"{b:>8s}" for b in BOARDS) + f"{'Overall':>9s}  10th-90th, Overall",
           row("A  each board at the better of its two uploads", r["A"], r["A"]["band"]["Overall"]),
           row("   A as of 2026-11-04 (c), before K2a", r["A_as_of_c"], r["A_as_of_c"]["band"]["Overall"]),
           row("B  first round: T1/heart/embryo s1, T3's s2", r["B"], r["B"]["band"]["Overall"]),
           row("   B as of 2026-11-04 (c), before K2a", r["B_as_of_c"], r["B_as_of_c"]["band"]["Overall"]),
           row("outlook of 2026-11-04 (c), docs/PLAN.md:990", res["outlook"], res["outlook"]["range"]),
           row("A as of (c) minus the outlook", res["A_as_of_c_minus_outlook"], signed=True),
           row("A minus the outlook", res["A_minus_outlook"], signed=True), ""]
    for name in ("A", "B"):
        out += _wrap(f"10th-90th by board, {name}: ",
                     [f"{b} [{r[name]['band'][b][0]:.1f}, {r[name]['band'][b][1]:.1f}]" for b in BOARDS])
    st = res["steps"]
    t1, means = st["T1"], st["T3_stand_in_means"]
    out += [f"steps: T1 is the grid under the {st['prior']} prior (A {t1['grid_A']:.2f}, B {t1['grid_B']:.2f}) plus "
            f"M1's gain over B0 on its backtest ({t1['m1_minus_b0']:+.2f})",
            "       heart moved by b0 on its backtest minus hold on the released truth: "
            + ", ".join(f"{h['task'].split('/')[1]} {h['offset']:+.2f}" for h in st["heart"]),
            "       T3 stand-in means: " + ", ".join(f"{k} {v:.2f}" for k, v in means.items()),
            "       T3 s1 weights as of (c): " + ", ".join(f"{k} {v}" for k, v in st["T3_weights_as_of_c"].items())
            + "; with K2a: " + ", ".join(f"{k} {v}" for k, v in st["T3_weights_with_k2a"].items())]
    for b, d in res["sensitivity_to_stated_priors"].items():
        out += _wrap(f"stated priors, {b} A/B: ", [f"{n} {v['A']:.2f}/{v['B']:.2f}" for n, v in d.items()])
    diffs = ", ".join(f"{b} {res['A_as_of_c_minus_outlook'][b]:+.2f}" for b in BOARDS)
    out.append(f"A as of (c) reproduces the outlook within {res['agree_within']} point per board: "
               f"{'yes' if res['reproduces_outlook'] else 'no'} ({diffs})")
    g = res["gate"]
    yn = {True: "yes", False: "no"}
    out += [f"the captain's gate threshold ({GATE}): {g['threshold']}; above it: A {yn[g['A_above']]}, "
            f"A as of (c) {yn[g['A_as_of_c_above']]}, B {yn[g['B_above']]}", f"  ({g['note']})",
            f"{len(res['checks'])} self-checks passed; {len(res['inputs_sha256'])} committed files read"]
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", metavar="PATH", help="also write the full result, checks and input hashes as JSON")
    a = ap.parse_args(argv)
    res = estimate()
    print(table(res))
    if a.json:
        out = Path(a.json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
