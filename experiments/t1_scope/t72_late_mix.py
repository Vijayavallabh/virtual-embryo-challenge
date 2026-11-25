#!/usr/bin/env python
"""Row T72 (experiments/harness/TRIES.md; a labelled diagnostic, twelve scenarios, one counted arm, one truth: half of
E10.5's library 3, zero time gap; never an adoption result): the practice read that binds row R54's late-mix file
T1_s2rcMSB_LM16. It prices only the payoff's shape: where the late-mix draw's lead over library 3's own mix turns
positive, and whether it survives where the blood premise fails. No released task holds E12.5's mix.

Geometry, truth half, reference, floor, protocol and seeds 0-9: T30's (experiments/t1_scope/t30_pricing.py, which
imports T27's kind_gap.py), imported unchanged, depth 1.0 only (no thinning anywhere).

Scenarios (twelve): the truth half re-drawn by quota at a cardiomyocyte share of heart cells in {own (the half's own,
about 0.24), 0.30, 0.33, LM16's, 0.45, 0.60} x a blood share of chamber cells in {own, LM16's}, LM16's two shares
being row N64's (shares.json: cm_share_of_heart, blood_share_of_chamber):
  blood own   t30.cm_truth_rows unchanged (cm own: the half itself, T30's truth at depth 1.0)
  blood LM16  the same quota draw at with_share(with_share(own mix, "cm", cm), "blood", LM16's blood,
              within=chamber) (cardiomyocytes first, then blood within the chamber, the rest in proportion), the
              non-chamber cells cut by the same factor, on its own stream [TRUTH_SEED, round(1000 cm), round(1000 b)]

Arms (n cells by the count rule; real cells of library 3's other half, half B; none from the truth half):
  lib3_base  T27's draw (kind_gap.arm_rows "lib3_base"), rerun, not counted: it must reproduce T30's stored per-seed
             scores at (own, own) within 0.02
  lib3_lm16  counted (one try): the same draw rule and the same stream ([seed, 3], lib3_base's) at half B's own mix
             moved to row N64's two shares by the builder's --late-mix rule (with_share: cardiomyocytes of heart
             cells, then blood of chamber cells, the other categories in proportion); quotas at the count, never
             weights. With row N64's shares at half B's own the two arms are one draw, byte for byte.

Bars, fixed in row T72 (d = lib3_lm16 - lib3_base per seed):
  (a) at the scenario carrying LM16's two shares, the 10-seed mean of d at least 2.0 with at least 8 of 10 seeds
      above 0;
  (b) the same at LM16's cardiomyocyte share with blood at own;
  (c) the break-even cardiomyocyte share at own blood (where the 10-seed mean of d turns positive for good, linearly
      between adjacent grid shares, each scenario's share as its truth realises it) at most the logit midpoint of
      library 3's own share and LM16's (row N64's shares.json) plus 0.03.
PASS only if (a), (b) and (c) all hold; otherwise STOP; an arm byte-equal to lib3_base in every seed of every scenario
reads STOP. VOID (nothing follows; the try counts): a scenario, seed or arm missing; results of more than one commit,
not at the row's commit, or from a dirty tree; lib3_base more than 0.02 from T30's stored per-seed score at (own,
own); row N64's shares more than 1e-9 from an independent recomputation (experiments/t1_scope/late_mix.py's second
path, on the same E10.5 file and manifest) or a shares file other than the one rows N64 and T72 name; an arm drawing
a truth-half cell. Reported, deciding nothing: every scenario's mean d and seeds ahead, the arm's repeat share, its
drawn median library over lib3_base's ((f)(3): a ratio of released cells, no count).

  run   --cm own|0.30|0.33|lm16|0.45|0.60 --blood own|lm16 --commit C   one scenario (the launcher starts the twelve
        together). Refused before anything loads unless row T72 is in TRIES.md naming this script, HEAD is C, the
        tracked tree is clean, the scenario file does not exist, and the committed shares file has the sha256 rows
        N64 and T72 both name. Prints no score.
  read  --commit C   the bars, after every void condition; writes t72.json and t72.txt; prints the outcome words and
        what follows only.
  smoke              before row T72 is written (no row needed, no try, no score, nothing written): the twelve
        scenarios' construction on the E10.5 file, as T30's was smoke-tested: per scenario the truth's cells and its two
        realised shares, and on seed 0's draw each arm's two realised shares, the arm's repeat share, its median library
        over lib3_base's and whether any arm draws a truth-half cell. Needs row N64's committed shares file (row N64
        naming exactly one `shares.json sha256 <hex>`, the file's).

Usage: through experiments/t1_scope/t72_late_mix.sh (one step per call, from the root of a checkout holding
data/raw/T1/E10.5_RNA.h5ad); the smoke test directly, `uv run --frozen --no-sync python -I -B
experiments/t1_scope/t72_late_mix.py smoke`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import re
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import kind_gap as kg  # noqa: E402  (T27's geometry and draws, unchanged)
import late_mix as lm  # noqa: E402  (row N64's step: its independent recomputation)
import t30_pricing as t30  # noqa: E402  (T30's truth re-draw, unchanged)

from vec import harness, t1_scope  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import RANKED_METRICS, load_boards  # noqa: E402
from vec.core import quotas  # noqa: E402
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample  # noqa: E402
from vec.prep import load_stage  # noqa: E402

ROW = "T72"
SCRIPT_PATH = "experiments/t1_scope/t72_late_mix.py"
TRIES = ROOT / "experiments" / "harness" / "TRIES.md"
RESULTS = HERE / "results" / "t72_late_mix"
SCENARIOS = RESULTS / "scenarios"
OUT_FILES = ("t72.json", "t72.txt")
T30_STORED = HERE / "results" / "t30_pricing" / "scenarios" / "t30_d1.0_cmown.json"
SHARES = lm.OUT / "shares.json"  # row N64's output
SHA_NAMED = re.compile(r"`shares\.json sha256 ([0-9a-f]{64})`")
CMS = ("own", "0.30", "0.33", "lm16", "0.45", "0.60")
BLOODS = ("own", "lm16")
BASE, ARM = "lib3_base", "lib3_lm16"
ARMS = (BASE, ARM)
TRUTH_SEED = 20261125  # the truth half's quota re-draw where blood moves (T30's own stream where it does not)
REPRO_TOL = 0.02  # points per seed, the tolerance of the ledger's reproduction gates
RECOMPUTE_TOL = lm.TOL  # 1e-9
BAR_LEAD = 2.0
MIN_AHEAD = 8
TRIES_COUNT = 1  # lib3_lm16; lib3_base is a rerun, not counted
LABEL = ("try T72, a labelled diagnostic on one truth (half of E10.5's library 3) at zero time gap, never an adoption "
         "result: row R54's late-mix draw (library 3's other half at row N64's shares) against library 3's own mix "
         "over twelve truth mixes")
FOLLOWS = {"PASS": ("What follows: a priced question to the captain with both readings (the level reading of row "
                    "N64's shares and the change transfer it reports); only then is T1_s2rcMSB_LM16 built, twice with "
                    "one sha256, by experiments/builds/t1_lm_m20_line.sh. Nothing is built or changed by this read."),
           "STOP": ("What follows: the late-mix line ends. No question goes to the captain, nothing is built, and "
                    "nothing is rerun under another row."),
           "VOID": ("What follows: nothing; the try counts, no file is built and nothing is rerun under this row.")}
_D: dict = {}


# ---------------------------------------------------------------------------------------------------------------------
# Row N64's shares


def named_shas(row: str, first_cell_only: bool, tries: Path | None = None) -> set[str] | None:
    """The distinct `shares.json sha256 <hex>` strings in the one row `row` (its first cell only, or every cell);
    None unless the ledger holds exactly one such row."""
    tries = TRIES if tries is None else tries
    lines = [ln for ln in tries.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    if len(lines) != 1:
        return None
    text = lines[0][2:].split("|", 1)[0] if first_cell_only else lines[0]
    return set(SHA_NAMED.findall(text))


def shares_sha_named(tries: Path | None = None) -> str | None:
    """The one sha256 both rows name: row N64 (any cell) and row T72 (its first cell), each exactly once; else None."""
    n64, t72 = named_shas(lm.ROW, False, tries), named_shas(ROW, True, tries)
    if not n64 or not t72 or len(n64) != 1 or n64 != t72:
        return None
    return next(iter(n64))


def committed_bytes(path: Path) -> bytes | None:
    rel = Path(path).resolve().relative_to(ROOT).as_posix()
    out = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT, capture_output=True)
    return out.stdout if out.returncode == 0 else None


def load_shares_file(sha: str, path: Path | None = None, committed=None) -> dict:
    """Row N64's shares file, committed at HEAD and the same bytes in the tree, with the sha256 the rows name."""
    path = SHARES if path is None else Path(path)
    raw = (committed_bytes if committed is None else committed)(path)
    if raw is None or not path.is_file() or path.read_bytes() != raw:
        sys.exit(f"{ROW} is refused: row {lm.ROW}'s shares file is not committed at HEAD as it is in the tree")
    if hashlib.sha256(raw).hexdigest() != sha:
        sys.exit(f"{ROW} is refused: the shares file's sha256 is not the one rows {lm.ROW} and {ROW} name")
    rec = json.loads(raw)
    if rec.get("row") != lm.ROW or rec.get("late_mix") != lm.LATE_MIX:
        sys.exit(f"{ROW} is refused: the shares file is not row {lm.ROW}'s output for {lm.LATE_MIX}")
    return rec


# ---------------------------------------------------------------------------------------------------------------------
# Truths and arms


def scenario_values(cm: str, blood: str, shares: dict) -> tuple[float | None, float | None]:
    """A scenario's cardiomyocyte share of heart cells and blood share of chamber cells (None: the truth half's own)."""
    if cm not in CMS or blood not in BLOODS:
        raise ValueError(f"--cm must be one of {CMS} and --blood one of {BLOODS}")
    c = None if cm == "own" else float(shares["cm_share_of_heart"]) if cm == "lm16" else float(cm)
    b = None if blood == "own" else float(shares["blood_share_of_chamber"])
    return c, b


def truth_rows(d: dict, cm: float | None, blood: float | None) -> np.ndarray:
    """The truth half re-drawn at a cardiomyocyte share of heart cells and a blood share of chamber cells: T30's own
    rows where blood is the half's own (`t30.cm_truth_rows`; the half itself where both are), else the same quota draw
    at the moved mix (module docstring) on its own stream."""
    if blood is None:
        return d["half_a"] if cm is None else t30.cm_truth_rows(d, cm)
    cat, half, chamber = d["cat"], d["half_a"], d["chamber"]
    counts = {c: int((cat[half] == c).sum()) for c in chamber}
    total = sum(counts.values())
    own = {c: counts[c] / total for c in chamber}
    target = own if cm is None else t1_scope.with_share(own, "cm", cm)
    target = t1_scope.with_share(target, "blood", blood, within=list(target))
    n_ch = int(min(counts[c] / target[c] for c in chamber if target[c] > 0))
    rng = np.random.default_rng([TRUTH_SEED, 0 if cm is None else int(round(cm * 1000)), int(round(blood * 1000))])
    keep = []
    for c in chamber:
        pool = half[cat[half] == c]
        keep.append(pool[rng.choice(len(pool), min(len(pool), int(n_ch * target[c])), replace=False)])
    rest = half[~np.isin(cat[half], chamber)]
    keep.append(rest[rng.choice(len(rest), int(round(len(rest) * n_ch / total)), replace=False)])
    return np.sort(np.concatenate(keep))


def lm16_mix(lib3_mix: dict, shares: dict) -> dict:
    """Half B's own mix moved by the builder's --late-mix rule to row N64's two shares."""
    m = t1_scope.with_share(dict(lib3_mix), "cm", float(shares["cm_share_of_heart"]))
    return t1_scope.with_share(m, "blood", float(shares["blood_share_of_chamber"]), within=list(m))


def arm_rows(arm: str, n: int, seed: int) -> np.ndarray:
    """lib3_base: kind_gap's draw. lib3_lm16: kind_gap's lib3_base draw rule and stream at the moved mix."""
    if arm == BASE:
        return kg.arm_rows(BASE, n, seed)
    d = kg._D
    rng = np.random.default_rng([seed, kg.ARMS.index(BASE)])
    q = quotas(d["lm16_mix"], n)
    return np.sort(np.concatenate([kg._draw(d["pools"]["b"][c], int(q.get(c, 0)), rng) for c in d["chamber"]]))


def shares_of(d: dict, rows: np.ndarray) -> tuple[float, float]:
    """A set of rows' cardiomyocyte share of heart cells and blood share of chamber cells."""
    heart = list(t1_scope.load_spec().heart)
    tc = {c: int((d["cat"][rows] == c).sum()) for c in d["chamber"]}
    h = sum(tc[c] for c in heart)
    return tc["cm"] / max(h, 1), tc.get("blood", 0) / max(h + tc.get("blood", 0), 1)


def prepare(data_dir: Path, cm: str, blood: str, shares: dict) -> dict:
    return scenario(kg.prepare(data_dir, 1.0), cm, blood, shares)  # T27's split, pools, mixes, the truth as released


def scenario(d0: dict, cm: str, blood: str, shares: dict) -> dict:
    """One scenario on kind_gap's prepared geometry `d0` (left as it is): the truth's rows and board, the arm's mix."""
    c, b = scenario_values(cm, blood, shares)
    d = dict(d0)
    rows = truth_rows(d, c, b)
    boards = load_boards()
    d["board"] = replace(boards[kg.TWIN], key=f"pseudo:T1/E10.5 late-mix pricing cm{cm} blood{blood}", split="pseudo",
                         anchors={}, **harness.pseudo_counts(len(rows), len(d["ref_rows"])))
    d.update(x_truth=d["X"][rows], truth_rows=rows, cm=cm, blood=blood, lm16_mix=lm16_mix(d["lib3_mix"], shares))
    d["truth_cm_of_heart"], d["truth_blood_of_chamber"] = shares_of(d, rows)
    d["truth_cells"] = int(len(rows))
    return d


def score_seed(seed: int) -> dict:
    d = kg._D
    t0 = time.time()
    X, board = d["X"], d["board"]
    target, ref = Stage(X=d["x_truth"]), Stage(X=X[d["ref_rows"]])
    f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, seed), target, ref, seed=seed, anchors="computed",
                    build=harness.BUILD, protocol=harness.PROTOCOL)
    if abs(f.score - 50.0) > 1e-9 or set(RANKED_METRICS["T1"]) - set(f.anchors):
        raise ValueError(f"seed {seed}: the floor does not score 50")
    n = harness.count_for(board, len(d["ref_rows"]))
    row = {"seed": seed, "scores": {}, "skills": {}, "median_library": {}, "repeat_share": {}, "rows_sha256": {},
           "draws_truth_cell": {}}
    for arm in ARMS:
        rows = arm_rows(arm, n, harness.BUILD_SEED_OFFSET + seed)
        r = score_board(board, Stage(X=X[rows]), target, ref, seed=seed, anchors=f.anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        row["scores"][arm], row["skills"][arm] = r.score, r.skills
        row["median_library"][arm] = float(np.median(d["tot"][rows]))
        row["repeat_share"][arm] = float(1.0 - len(np.unique(rows)) / len(rows))
        row["rows_sha256"][arm] = hashlib.sha256(np.ascontiguousarray(rows, dtype=np.int64).tobytes()).hexdigest()
        row["draws_truth_cell"][arm] = bool(np.isin(rows, d["half_a"]).any())
    row["seconds"] = round(time.time() - t0, 1)
    print(f"{d['cm']} {d['blood']} seed {seed}: {row['seconds']:.0f} s", flush=True)  # no score
    return row


def scenario_name(cm: str, blood: str) -> str:
    return f"t72_cm{cm}_blood{blood}.json"


def smoke(d0: dict, shares: dict) -> list[dict]:
    """The twelve scenarios' construction without scoring (module docstring): per scenario the truth's cells and two
    realised shares, and on seed 0's draw each arm's two realised shares, the arm's repeat share, its median library
    over lib3_base's and whether an arm draws a truth-half cell."""
    out = []
    for cm in CMS:
        for blood in BLOODS:
            d = scenario(d0, cm, blood, shares)
            kg._D.clear()
            kg._D.update(d)
            n = harness.count_for(d["board"], len(d["ref_rows"]))
            drawn = {arm: arm_rows(arm, n, harness.BUILD_SEED_OFFSET) for arm in ARMS}
            arm_cm, arm_blood = shares_of(d, drawn[ARM])
            base_cm, base_blood = shares_of(d, drawn[BASE])
            out.append({"cm": cm, "blood": blood, "truth_cells": d["truth_cells"],
                        "truth_cm_of_heart": d["truth_cm_of_heart"],
                        "truth_blood_of_chamber": d["truth_blood_of_chamber"],
                        "base_cm_of_heart": base_cm, "base_blood_of_chamber": base_blood, "lm16_cm_of_heart": arm_cm,
                        "lm16_blood_of_chamber": arm_blood,
                        "repeat_share": float(1.0 - len(np.unique(drawn[ARM])) / len(drawn[ARM])),
                        "median_library_over_base": float(np.median(d["tot"][drawn[ARM]])
                                                          / np.median(d["tot"][drawn[BASE]])),
                        "draws_truth_cell": bool(any(np.isin(drawn[a], d["half_a"]).any() for a in ARMS))})
    return out


def cmd_smoke(a) -> None:
    named = named_shas(lm.ROW, False)
    if not named or len(named) != 1:
        sys.exit(f"{ROW}'s smoke test is refused: row {lm.ROW} does not name exactly one `shares.json sha256 <hex>`")
    shares = load_shares_file(next(iter(named)))
    for r in smoke(kg.prepare(Path(a.data), 1.0), shares):
        print(f"cm {r['cm']}, blood {r['blood']}: truth {r['truth_cells']} cells, cardiomyocytes "
              f"{r['truth_cm_of_heart']:.3f} of heart, blood {r['truth_blood_of_chamber']:.3f} of chamber; seed 0's "
              f"draws: {BASE} {r['base_cm_of_heart']:.3f} / {r['base_blood_of_chamber']:.3f}, {ARM} "
              f"{r['lm16_cm_of_heart']:.3f} / {r['lm16_blood_of_chamber']:.3f}, repeats {r['repeat_share']:.3f}, "
              f"median library over {BASE}'s {r['median_library_over_base']:.4f}, a truth-half cell drawn "
              f"{r['draws_truth_cell']}", flush=True)


def refuse_run(commit: str | None, path: Path) -> str:
    """The refusals before anything loads; returns the shares sha256 the rows name."""
    first = lm.row_first_cell(ROW, TRIES)
    if first is None or not first.startswith(f"{ROW} (") or not lm.names_path(first, SCRIPT_PATH):
        sys.exit(f"{ROW} is refused: experiments/harness/TRIES.md has no single row {ROW} whose first cell starts "
                 f"'{ROW} (' and names {SCRIPT_PATH}")
    if not commit or harness._git("rev-parse", "HEAD") != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the commit given (the full hash of row {ROW}'s commit)")
    if harness._git("status", "--porcelain", "--untracked-files=no") != "":
        sys.exit(f"{ROW} is refused: the tracked tree is not clean")
    if path.exists():
        sys.exit(f"{ROW} is refused: {path.name} exists (a scenario is never run twice)")
    sha = shares_sha_named()
    if sha is None:
        sys.exit(f"{ROW} is refused: rows {lm.ROW} and {ROW} do not each name the one same `shares.json sha256 <hex>`")
    return sha


def cmd_run(a) -> None:
    path = SCENARIOS / scenario_name(a.cm, a.blood)
    sha = refuse_run(a.commit, path)
    shares = load_shares_file(sha)
    ref = load_stage(Path(a.data) / lm.E105, lineage=False)  # the independent recomputation of row N64's shares
    late = json.loads(lm.LATE_MANIFEST.read_bytes())
    worst = lm.largest_difference(shares, lm.recompute(lm.own_shares_direct(ref), late["label_counts"]))
    del ref
    kg._D.update(prepare(Path(a.data), a.cm, a.blood, shares))
    d = kg._D
    if a.jobs > 1:
        with mp.get_context("fork").Pool(a.jobs) as pool:
            rows = pool.map(score_seed, list(SEEDS))
    else:
        rows = [score_seed(s) for s in SEEDS]
    res = {"try": ROW, "label": LABEL, "commit": a.commit, "tree_clean": True, "cm": a.cm, "blood": a.blood,
           "truth_cells": d["truth_cells"], "truth_cm_of_heart": d["truth_cm_of_heart"],
           "truth_blood_of_chamber": d["truth_blood_of_chamber"], "seeds": list(SEEDS), "arms": list(ARMS),
           "protocol": harness.PROTOCOL, "build": harness.BUILD, "frac": SUBSAMPLE_FRAC,
           "library3_chamber_mix": d["lib3_mix"], "lm16_mix": d["lm16_mix"], "shares_sha256": sha,
           "recomputation": {"largest_difference": worst, "tolerance": RECOMPUTE_TOL},
           "rows": sorted(rows, key=lambda r: r["seed"])}
    SCENARIOS.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1) + "\n")
    print(f"wrote {path.name}", flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def logit(p: float) -> float:
    return math.log(p / (1.0 - p))


def expit(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def break_even(points) -> float:
    """points: (share, mean d) at own blood, any order. The share where d turns positive for good, linearly between
    adjacent grid shares: above the last grid share where d <= 0, the crossing to the next; the lowest grid share if
    d > 0 everywhere; +inf if d <= 0 at the largest share."""
    pts = sorted((float(x), float(y)) for x, y in points)
    if pts[-1][1] <= 0:
        return math.inf
    last = max((i for i, (_, y) in enumerate(pts) if y <= 0), default=None)
    if last is None:
        return pts[0][0]
    (x0, y0), (x1, y1) = pts[last], pts[last + 1]
    return x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0)


def diff(r: dict) -> dict:
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    d = np.array([x["scores"][ARM] - x["scores"][BASE] for x in rows], dtype=np.float64)
    ratio = [x["median_library"][ARM] / x["median_library"][BASE] for x in rows]
    return {"mean": float(d.mean()), "sd": float(d.std(ddof=1)) if d.size > 1 else 0.0,
            "ahead": int((d > 0).sum()), "behind": int((d < 0).sum()), "of": int(d.size),
            "range": [float(d.min()), float(d.max())],
            "mean_scores": {a: float(np.mean([x["scores"][a] for x in rows])) for a in ARMS},
            "repeat_share": float(np.mean([x["repeat_share"][ARM] for x in rows])),
            "median_library_over_base": float(np.mean(ratio)),
            "truth_cm_of_heart": float(r["truth_cm_of_heart"]),
            "truth_blood_of_chamber": float(r["truth_blood_of_chamber"])}


def guards(loaded: dict, stored: dict | None, sha: str | None, row_commit: str | None) -> list[str]:
    """Every void condition of row T72 (module docstring); empty when none fired."""
    why = []
    for cm in CMS:
        for blood in BLOODS:
            r = loaded.get((cm, blood))
            if r is None:
                why.append(f"scenario cm {cm}, blood {blood} is missing")
                continue
            rows = r.get("rows") or []
            if sorted(int(x["seed"]) for x in rows) != list(SEEDS) or any(
                    a not in x.get("scores", {}) or a not in x.get("rows_sha256", {}) for x in rows for a in ARMS):
                why.append(f"cm {cm}, blood {blood}: a seed or an arm is missing")
            if any(any((x.get("draws_truth_cell") or {}).get(a, True) for a in ARMS) for x in rows):
                why.append(f"cm {cm}, blood {blood}: an arm draws a truth-half cell")
            rec = (r.get("recomputation") or {}).get("largest_difference")
            if rec is None or not rec <= RECOMPUTE_TOL:
                why.append(f"cm {cm}, blood {blood}: row {lm.ROW}'s shares differ from the independent "
                           f"recomputation by {rec} (tolerance {RECOMPUTE_TOL:g})")
            if sha is None or r.get("shares_sha256") != sha:
                why.append(f"cm {cm}, blood {blood}: it ran on another shares file than the one the rows name")
    present = [r for r in loaded.values() if r is not None]
    if len({r.get("commit") for r in present}) > 1:
        why.append("the results come from more than one commit")
    if row_commit is not None and any(r.get("commit") != row_commit for r in present):
        why.append("a result did not run at the row's commit")
    if not all(r.get("tree_clean") is True for r in present):
        why.append("a result comes from a dirty tree")
    base = loaded.get(("own", "own"))
    if not stored or not stored.get("rows"):
        why.append("T30's stored scores at depth 1.0, cm own are missing")
    elif base is not None and not why:
        old = {int(x["seed"]): x["scores"][BASE] for x in stored["rows"]}
        worst = max((abs(x["scores"][BASE] - old[int(x["seed"])]) if int(x["seed"]) in old else math.inf
                     for x in base["rows"]), default=math.inf)
        if not worst <= REPRO_TOL:
            why.append(f"lib3_base does not reproduce T30's stored per-seed scores at (own, own) (largest difference "
                       f"{worst:.4g}, tolerance {REPRO_TOL:g})")
    return why


def idle(loaded: dict) -> bool:
    """lib3_lm16 is lib3_base's draw in every seed of every scenario (by the rows' sha256)."""
    return all(x["rows_sha256"][ARM] == x["rows_sha256"][BASE] for r in loaded.values() for x in r["rows"])


def t72_reading(loaded: dict, stored: dict | None, shares: dict | None, sha: str | None,
                row_commit: str | None = None) -> dict:
    """Row T72's reading: the void conditions first, then the bars."""
    why = guards(loaded, stored, sha, row_commit)
    if shares is None:
        why.append(f"row {lm.ROW}'s shares file was not read")
    if why:
        return {"outcome": "VOID", "reasons": why, "follows": FOLLOWS["VOID"]}
    by = {f"cm {cm}, blood {blood}": diff(loaded[(cm, blood)]) for cm in CMS for blood in BLOODS}
    a, b = by["cm lm16, blood lm16"], by["cm lm16, blood own"]
    own_cm, lm_cm = float(shares["own_cm_share_of_heart"]), float(shares["cm_share_of_heart"])
    mid = expit((logit(own_cm) + logit(lm_cm)) / 2.0)
    threshold = mid + lm.BAR_C_MARGIN
    points = [(by[f"cm {cm}, blood own"]["truth_cm_of_heart"], by[f"cm {cm}, blood own"]["mean"]) for cm in CMS]
    even = break_even(points)
    bars = {"a": bool(a["mean"] >= BAR_LEAD and a["ahead"] >= MIN_AHEAD),
            "b": bool(b["mean"] >= BAR_LEAD and b["ahead"] >= MIN_AHEAD), "c": bool(even <= threshold)}
    failed = []
    if not bars["a"]:
        failed.append(f"(a) at LM16's two shares the lead is not {BAR_LEAD:g} or more with {MIN_AHEAD} of 10 seeds "
                      "ahead")
    if not bars["b"]:
        failed.append(f"(b) at LM16's cardiomyocyte share with blood at own the lead is not {BAR_LEAD:g} or more with "
                      f"{MIN_AHEAD} of 10 seeds ahead")
    if not bars["c"]:
        failed.append("(c) the break-even cardiomyocyte share at own blood is above the logit midpoint plus "
                      f"{lm.BAR_C_MARGIN:g}")
    still = idle(loaded)
    if still:
        failed.append("lib3_lm16 is lib3_base's draw in every seed of every scenario")
    outcome = "PASS" if not failed else "STOP"
    return {"outcome": outcome, "failed": failed, "bars": bars, "by_scenario": by,
            "break_even": {"share": even, "logit_midpoint": mid, "threshold": threshold, "points": points,
                           "own_cm_share_of_heart": own_cm, "lm16_cm_share_of_heart": lm_cm},
            "idle": still, "follows": FOLLOWS[outcome]}


def outcome_lines(r: dict) -> list[str]:
    if r["outcome"] == "VOID":
        return [f"Outcome: {ROW} VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]]
    first = f"Outcome: {ROW} {r['outcome']}" + ("" if r["outcome"] == "PASS" else " (" + "; ".join(r["failed"]) + ")")
    return [first, r["follows"]]


def t72_text(res: dict) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        lines.append(f"d = {ARM} - {BASE} per seed (bars: (a) {r['bars']['a']}, (b) {r['bars']['b']}, (c) "
                     f"{r['bars']['c']}):")
        for k, p in r["by_scenario"].items():
            lines.append(f"  {k}: {p['mean']:+.2f} (sd {p['sd']:.2f}, ahead {p['ahead']} behind {p['behind']} of "
                         f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}); truth cardiomyocytes "
                         f"{p['truth_cm_of_heart']:.3f} of heart, blood {p['truth_blood_of_chamber']:.3f} of chamber; "
                         f"mean scores {BASE} {p['mean_scores'][BASE]:.2f}, {ARM} {p['mean_scores'][ARM]:.2f}; "
                         f"repeats {p['repeat_share']:.3f}; median library over {BASE}'s "
                         f"{p['median_library_over_base']:.4f}")
        e = r["break_even"]
        lines.append(f"break-even cardiomyocyte share at own blood {e['share']:.4f}; the logit midpoint of "
                     f"{e['own_cm_share_of_heart']:.4f} and {e['lm16_cm_share_of_heart']:.4f} is "
                     f"{e['logit_midpoint']:.4f}, the bar {e['threshold']:.4f}")
    lines.append(f"Tries counted: {res['tries']} (the arm {ARM}; {BASE} is a rerun)")
    lines.append(f"Commit {res['commit']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    first = lm.row_first_cell(ROW, TRIES)
    if first is None or not lm.names_path(first, SCRIPT_PATH):
        sys.exit(f"{ROW} is refused: no single row {ROW} names {SCRIPT_PATH}")
    if not a.commit or harness._git("rev-parse", "HEAD") != a.commit:
        sys.exit(f"{ROW} is refused: HEAD is not the row's commit given")
    if harness._git("status", "--porcelain", "--untracked-files=no") != "":
        sys.exit(f"{ROW} is refused: the tracked tree is not clean")
    there = [f for f in OUT_FILES if (RESULTS / f).exists()]
    if there:
        sys.exit(f"{ROW} is refused: {there} exist (a reading is never overwritten)")
    sha = shares_sha_named()
    shares = load_shares_file(sha) if sha is not None else None
    loaded = {}
    for cm in CMS:
        for blood in BLOODS:
            p = SCENARIOS / scenario_name(cm, blood)
            loaded[(cm, blood)] = json.loads(p.read_text()) if p.exists() else None
    stored = json.loads(T30_STORED.read_text()) if T30_STORED.exists() else None
    reading = t72_reading(loaded, stored, shares, sha, a.commit)
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT, "commit": a.commit,
           "shares_sha256": sha}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / OUT_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (RESULTS / OUT_FILES[1]).write_text(t72_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {RESULTS / OUT_FILES[1]})")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--cm", required=True, choices=CMS)
    r.add_argument("--blood", required=True, choices=BLOODS)
    r.add_argument("--commit", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--jobs", type=int, default=2)
    d = sub.add_parser("read")
    d.add_argument("--commit", required=True)
    s = sub.add_parser("smoke")
    s.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read, "smoke": cmd_smoke}[a.cmd](a)


if __name__ == "__main__":
    main()
