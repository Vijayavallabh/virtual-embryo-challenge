#!/usr/bin/env python
"""Try K30 (experiments/harness/TRIES.md; every step refuses without row K30 and ruling row R53, on an unclean tracked
tree or at a HEAD other than the commit given): the read to which the captain's named exception (ruling row R53, given
before row K30) binds one T3 upload: the frozen per-label depth step on C2-like cells (proposal P1 of the design round
of 2026-11-25, `dD`, as its cold challenge fixed it). A labelled analogue (rule 9, J1a) and a diagnostic (J5), not the
adoption rule. On the two Gata4 replicates, in row K18's four cells, it asks what the depth step gains where the
knockout sits at its run's depth (the plain cells: the truth at its run's depth, the reference and the pool at WT
E8.75's) and what it costs where it sits at its reference's (the dm cells: truth, reference and pool all at the
replicate's depth). It names at most one of three depth forms, or none.

Built on k18_built_files.py, imported unchanged: row K18's four cells (plain g1, plain g2, dm g1, dm g2), its stream
(`t3_s2_form`), seeds, barred spec and factor, plain build, depth-matched row, stamp, task-file writer, folder rule,
masking and guard distance. For a run, three of k18's module names are replaced inside this process only (`k18_arms`):
ARMS, build_arm and fold_factors, restored on exit, as rows K21 to K24 and K29 do it. Also imported unchanged:
k21_factor_forms.py (the seed parser, K18's committed rows at their pins, _read_json, _gated), k22_outside_forms.py (the
frozen-file check, the stamp check) and k7_wt_revert.py (_paired, MIN_AHEAD).

The arms. The base, `barred_form`, is K18's own build_arm for its barred_form (bs.construct with the barred spec's
keywords through prediction_inputs: phi_M at exponent 1, the central composition at D7m, no depth or stage step) on
stable_rng(seed, t3_s2_form); it is also the guard. Each depth arm is the same call with one keyword set, construct's
`ratio`, which t3_decompose.build hands rescale_depth after the factor on the same generator, so an arm and the base
share their rows and factor draws and differ by the depth step alone:
  dx1   the OTHER replicate's per-label ratios from the committed gata4_composition.csv (g1 takes rep2_depth_ratio, g2
        takes rep1_depth_ratio), put on the world's vocabulary by the builder's own depth_ratios (the labels whose
        ratio is not 1.0, as submission_parameters lists the pooled ones; every other label stays at 1.0)
  dxh   the same ratios ** 0.5
  dxg   one scalar for every label: the geometric mean of the other replicate's ratios, weighted by WT E8.75's cells per
        label (the table's wt_cells) over all labels, 1.0 where the frozen run measured none
In the dm cells the step goes on K18's matched pool by U3d's path for s1's stand-in (orth_diag.score_row's
dm:xrep:F+depth: t3_decompose.build with the matched pool as its source and the other replicate's ratio as `ratio`),
with K18's base in the place of U3d's factor and composition; each dm row file records how far its dx1 vector lies from
the live context's other-replicate ratio, U3d's S1.ratio (reported). There is no x2 cell. The table is read at its
pinned sha256 (refused otherwise); every ratio must be finite and positive and every weight finite and non-negative
with a positive sum. The vectors are made on each world's own vocabulary when the fold is prepared, so the guard part,
which scores the base alone, already makes them: a label the vocabulary lacks stops the guard part (a crash,
INCOMPLETE) before any counted arm is scored. The task (g1 or g2) whose ratios an arm takes is set for each build by
this script's wrappers (the plain build from the harness's task key, the depth-matched row from its target's name) and
cleared after it; an arm built with no task set raises. Each build's counts are observed as build_submission's
observe_draws observes draw_rows (the frozen functions called unchanged): every task file records, per arm, the mean
counts per built cell and whether each depth arm's counts before the step are byte for byte the base's built counts.

S3's world weight is named in row K30's first cell, exactly one of `--s3-weight 0.5` (both worlds, half each) and
`--s3-weight plain` (the plain cells only, the captain's ruling that the board keeps its best upload), and ruling row
R53's first cell must name the same one ("row K30's world weight is `plain`"); run and read take it as --s3-weight and
refuse another, every task file records it and the read voids a file run under another.

Steps (each with --commit, the full hash of the run commit, which must be HEAD; outputs in the run folder R):
  run --part guard --world plain --task g1|g2 [--seeds 0-4] --s3-weight W --out R   the guard part: the base alone,
  run --part guard --world dm [--workers 8 --threads 2] --s3-weight W --out R       into R/guard/; refused with K18's
                          committed files off their pins
  guard --runs R          the guard check (R/guard.json, made once): the base within K17's tolerance (ten-seed mean 0.1,
                          every seed 0.5 board points) of K18's committed barred_form in all four cells. Exit 6 when it
                          fails (VOID with one counted try, the base; the run step never starts)
  run [--part score] --world plain|dm ... --s3-weight W --out R   a scoring part: the base and the three depth arms;
                          refused without a passing guard.json of this commit and machine. THE THREE TRIES COUNT FROM
                          THE START OF THE RUN STEP: the first scoring part to pass its refusals writes
                          R/run_started.json. Either part refuses (exit 4) once a reading is in
                          experiments/t3_decomp/results/k30_depth_read: no arm is scored again, whatever the run folder
  read --runs R --s3-weight W [--out DIR]   the fixed reading: DIR/reading.json and DIR/reading.txt (DIR defaults to
                          experiments/t3_decomp/results/k30_depth_read), never over an existing one; reading.txt opens
                          with the outcome in fixed words; read prints the outcome line and the decision table
run and guard refuse at or after the cutoff, 2026-11-25 13:00 UTC (exit 7): nothing is written then, so a folder left
INCOMPLETE at the cutoff stays so, and read, which runs at any time, writes its VOID reading.

Layout of R: guard/plain_<task>/seed_<s>.json, guard/dm/<task>_seed_<s>.json, guard.json, run_started.json,
plain_<task>/seed_<s>.json, dm/<task>_seed_<s>.json: one task file per seed and cell, each with its <file>.prov.json
(k18.write_task_file); dm/setup.json and guard/dm/setup.json. A rerun at the same commit on the same machine skips the
seeds already written; a folder holding a file of another commit or machine is refused (k18.check_folder). No log holds
a score: progress lines only, K18's timing line per dm row and the fixed self-check line; a crash prints its traceback
with every decimal number masked (k18.masked).

The reading (fixed in row K30; `reading`, `decide`):
  INCOMPLETE  guard.json, a guard file, or (after a passing guard) a scoring file missing: not a reading; read opens no
              seed file, writes nothing and exits 3. A read made at or after 2026-11-25 13:00 UTC with a file still
              missing writes instead the VOID reading "INCOMPLETE past its cutoff".
  VOID        checked first, and it wins (no STOP clause is applied): the guard failed at the guard step (one counted
              try, the base); a seed file that misses an arm (the plain ones the floor too) or is not its part's,
              task's or seed's; mixed commits or machines; a dirty tree; K18's committed task files off their pins; the
              guard step's base, read again, off K18's rows; the run step's base off the guard step's beyond K17's
              tolerance; a seed file whose depth record is not the pinned table's, or that was run under another S3
              weight. Outcome "Try K30: VOID: [reasons]; no file." The count: three (dx1, dxh, dxg) once the run step
              has started, void or not; the base counts one try only if it fails its guard (at the guard step: VOID
              with that one try alone; in the run step: four); a malformed or missing file adds none.
  Per cell, paired by seed: D_a = arm - barred_form. AHEAD: mean above +1.0 with at least 8 of 10 seeds positive (a
              seed at exactly 0 is not positive).
  EV_a        w x mean(D_a plain g1, plain g2) + (1 - w) x mean(D_a dm g1, dm g2), w 1/2 for `0.5` and 1 for `plain`
              (the dm cells then decide nothing), in exact fractions of the cells' float means.
  An arm passes S1 when it is AHEAD in both plain cells, and passes S2 when in neither plain cell it has more
              severity-gated seeds (r2_gate, beta_gate, no_response) than barred_form plus 1. The arms that may be
              named are those passing S1 and S2.
  STOP        any of: S1, no arm is AHEAD in both plain cells; S2, an arm passes S1 but none that passes S1 passes S2
              (the arm to be named has more severity-gated seeds than barred_form plus 1 in a plain cell); S3, the best
              EV_a among the arms that may be named is +1.0 or less; S4 (nothing changes), |mean D_a| is below 0.3 in
              both plain cells for every arm. Outcome "Try K30: STOP ([clauses]): no file; C2 stays."
  NAME        otherwise: the arm with the largest EV_a among those passing S1 and S2 (an exact tie goes to the first of
              dx1, dxh, dxg). Outcome "Try K30: NAME <arm>: the file may be built twice and uploaded under the
              captain's bound yes (R53)."
  REPORTED, deciding nothing: severity-gated seeds per arm and cell (all four cells), per-term skill differences, mean
              counts per built cell against barred_form's and the pre-step identity, each cell's arms (mean, sd,
              skills, severity states), dx1 - dxh, dx1 - dxg, the base against the floor (plain), gate flips against the
              base, the depth record (the table's sha256, the pairing, the global ratios, the distance to U3d's ratio).
              K29's gate before upload is not used.
Nothing changes after any score: the pairing, the exponents 1 and 0.5, the weights, +1.0, 8 of 10, the S3 weight, 0.3,
plus 1, the tie order, the seeds, 0.1, 0.5 and the cutoff of 2026-11-25 13:00 UTC.

Exit codes: 0 done; 1 a crash (masked traceback; seeds written so far stay); 2 refused (no row K30 or R53 as above, an
unclean tree, another HEAD, another S3 weight, the frozen T3 files off their pins, K18's committed files off their
pins, the table off its pin or not well formed, a scoring part without the passing guard.json it needs); 3 INCOMPLETE;
4 a reading or guard.json is there already (never overwritten; run: no arm is scored again); 6 the guard failed; 7 run
or guard refused at or after the cutoff.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
import tempfile
import time
import traceback
from contextlib import contextmanager
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (_paired, MIN_AHEAD, FULL_PANEL)
import k18_built_files as k18  # noqa: E402  (worlds, tasks, stream, barred spec and factor, builds, scoring; unchanged)
import k21_factor_forms as k21  # noqa: E402  (seed parser, K18's committed rows, _read_json, _gated; unchanged)
import k22_outside_forms as k22  # noqa: E402  (the frozen-file check, the stamp check; unchanged)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

SCRIPT = "experiments/t3_decomp/k30_depth_read.py"
LAUNCHER = "experiments/t3_decomp/run_k30.sh"
ROW_LABEL = "K30"
ROW_START = "K30 ("  # row K30's first cell starts so, names SCRIPT and LAUNCHER and one S3 weight
RULING = "R53"  # the captain's named exception, given before row K30
RULING_START = "R53 (the captain's "  # its first cell starts so
RULING_WEIGHT = "row K30's world weight is `{}`"  # and names row K30's world weight so
WEIGHT_OPTION = "--s3-weight"
S3_WEIGHTS = {"0.5": Fraction(1, 2), "plain": Fraction(1)}  # the plain cells' weight in EV; the dm cells take the rest
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"  # looked up at call time (row_present)
STREAM = k18.STREAM  # "t3_s2_form"
TASKS = k18.TASKS  # g1, g2
WORLDS = k18.WORLDS  # plain, dm
CELLS = ("plain g1", "plain g2", "dm g1", "dm g2")
PLAIN_CELLS = ("plain g1", "plain g2")
DM_CELLS = ("dm g1", "dm g2")

BASE = k18.BARRED  # "barred_form": the base and the guard, K18's own arm
DEPTH_ARMS = ("dx1", "dxh", "dxg")  # the three counted tries, in the tie order
ARMS = (BASE, *DEPTH_ARMS)
COUNTED = DEPTH_ARMS
EXPONENT = {"dx1": 1.0, "dxh": 0.5}  # the per-label forms; dxg is the global one
GLOBAL_ARM = "dxg"
FILE_OF = {"dx1": "T3_c2_s2_dd.h5ad", "dxh": "T3_c2_s2_ddh.h5ad", "dxg": "T3_c2_s2_dg.h5ad"}  # R53's names
REPORTED_PAIRS = (("dx1", "dxh"), ("dx1", "dxg"))

# The depth table and the pairing: each task takes the OTHER replicate's column
TABLE = HERE / "results" / "gata4" / "gata4_composition.csv"
TABLE_REL = "experiments/t3_decomp/results/gata4/gata4_composition.csv"
TABLE_SHA256 = "6d31f2f93cfc289c544ade5d3074aac9552ac772cc7b81666111bef1724bfcf1"
OTHER = {"g1": "rep2", "g2": "rep1"}
OTHER_TASK = {"g1": "g2", "g2": "g1"}  # whose truth the other replicate is
COLUMN = {n: f"{r}_depth_ratio" for n, r in OTHER.items()}
WEIGHTS = "wt_cells"  # WT E8.75's cells per label, the frozen run's vocabulary counts

# The reading's constants
NAME_GAIN = 1.0  # AHEAD: mean above +1.0 ...
MIN_AHEAD = k7.MIN_AHEAD  # ... with at least 8 of 10 seeds positive
MAX_EXTRA_GATED = 1  # S2: "barred_form plus 1"
EV_BAR = Fraction(1)  # S3: the best EV_a +1.0 or less stops
NULL_BAND = 0.3  # S4: |mean D_a| below 0.3 in both plain cells for every arm
GUARD_MEAN, GUARD_SEED = k18.GUARD_MEAN, k18.GUARD_SEED  # 0.1 and 0.5 (K17's tolerance)
DEPTH_REL = 1e-12  # the read's check of a recorded vector against the pinned table (an ulp of a power may differ)
CUTOFF = datetime(2026, 11, 25, 13, 0, tzinfo=UTC)
CUTOFF_TEXT = "2026-11-25 13:00 UTC"
K18_TASKS = k21.K18_TASKS  # experiments/t3_decomp/results/k18_built_files/tasks
K18_SHA256 = {"plain_g1.json": "5ec09ef583862ab47ab078e72c7aed807c89b4d9c234c57e3b0727ff9e460ed9",
              "plain_g2.json": "cfc5ede7e7185ba5e16d4774d975b9e140e60bbd29736fb85c0fdb8e99f4634f",
              "dm.json": "1ef3e26347c65f915c1630130da02a655c155b912392180001cba315fe135b60"}
RESULTS = HERE / "results" / "k30_depth_read"
GUARD_FILE, GUARD_DIR, RUN_MARK = "guard.json", "guard", "run_started.json"
EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_GUARD, EXIT_LATE = 2, 3, 4, 6, 7
VOID_LINE = "Try K30: VOID: {}; no file."
STOP_LINE = "Try K30: STOP ({}): no file; C2 stays."
NAME_LINE = "Try K30: NAME {arm}: the file may be built twice and uploaded under the captain's bound yes (R53)."
LATE = "INCOMPLETE past its cutoff"
CLAUSES = {"S1": "S1: no arm is AHEAD in both plain cells",
           "S2": "S2: no arm AHEAD in both plain cells keeps its severity-gated seeds within barred_form's plus 1 in "
                 "both plain cells",
           "S3": "S3: the best EV among the arms that may be named is +1.0 or less",
           "S4": "S4: nothing changes (|mean D| below 0.3 in both plain cells for every arm)"}
_K30: dict = {}  # "table"; "task" (the build in progress, set and cleared per build); "counts" (its builds' figures)


class Refused(Exception):
    """A refusal of the whole step (exit status 2)."""


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


def _now() -> datetime:
    return datetime.now(UTC)


def past_cutoff() -> bool:
    return _now() >= CUTOFF


# ---------------------------------------------------------------------------------------------------------------------
# The ledger rows and the preconditions


def _names(cell: str, path: str) -> bool:
    """`cell` names `path` as a whole path (not inside a longer one)."""
    return re.search(rf"(?<![\w./-]){re.escape(path)}(?![\w./-])", cell) is not None


def _first_cell(lines: list[str], rid: str) -> str:
    """The first cell of the one line opening with '| <rid> ' (refused unless there is exactly one)."""
    rows = [ln for ln in lines if ln.startswith(f"| {rid} ")]
    if len(rows) != 1:
        refuse(f"TRIES.md does not hold exactly one row {rid} ({len(rows)} lines open with '| {rid} '); the row is "
               "committed before any step")
    return rows[0].split("|")[1]


def row_present(path: Path | None = None) -> dict:
    """Rows K30 and R53, as the launcher's gate checks them (exit 2 otherwise): TRIES.md holds exactly one line opening
    with '| K30 ', whose first cell starts 'K30 (', names this script and its launcher as whole paths, and names exactly
    one S3 weight, `--s3-weight 0.5` or `--s3-weight plain`, as a whole backticked string (a backticked string holding
    --s3-weight with any other value refuses); and exactly one line opening with '| R53 ', whose first cell starts
    "R53 (the captain's " and names row K30's world weight as "row K30's world weight is `<weight>`", the same one.
    Returns {"id": "K30", "weight": <weight>, "ruling": "R53"}."""
    path = TRIES_MD if path is None else Path(path)
    try:
        lines = Path(path).read_text().splitlines()
    except OSError:
        lines = []
    first = _first_cell(lines, ROW_LABEL)
    if not first.startswith(f" {ROW_START}") or not all(_names(first, p) for p in (SCRIPT, LAUNCHER)):
        refuse(f"row {ROW_LABEL}'s first cell does not start '{ROW_START}' and name {SCRIPT} and {LAUNCHER}")
    named = {s for s in re.findall(r"`([^`]*)`", first) if WEIGHT_OPTION in s}
    weights = [w for w in S3_WEIGHTS if f"{WEIGHT_OPTION} {w}" in named]
    if len(named) != 1 or len(weights) != 1:
        refuse(f"row {ROW_LABEL}'s first cell does not name exactly one S3 weight, `{WEIGHT_OPTION} 0.5` or "
               f"`{WEIGHT_OPTION} plain`")
    weight = weights[0]
    ruling = _first_cell(lines, RULING)
    if not ruling.startswith(f" {RULING_START}"):
        refuse(f"row {RULING}'s first cell does not start \"{RULING_START}\"")
    said = [w for w in S3_WEIGHTS if RULING_WEIGHT.format(w) in ruling]
    if said != [weight]:
        refuse(f"row {RULING}'s first cell does not name row {ROW_LABEL}'s world weight as row {ROW_LABEL} does "
               f"({weight})")
    return {"id": ROW_LABEL, "weight": weight, "ruling": RULING}


DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing written; INCOMPLETE, not VOID")


def head_ok(st: dict, commit: str) -> None:
    if not commit or st["commit"] != commit:
        refuse(f"HEAD is {st['commit']}, not the run commit {commit}")


def weight_ok(row: dict, weight: str) -> None:
    if weight != row["weight"]:
        refuse(f"{WEIGHT_OPTION} {weight} is not the S3 weight row {ROW_LABEL} names ({row['weight']})")


LATE_REFUSAL = (f"past the cutoff ({CUTOFF_TEXT}): run and guard write nothing now, so a folder left INCOMPLETE stays "
                f"so; read writes the VOID reading '{LATE}'")


def preconditions(commit: str) -> tuple[dict, dict]:
    """Before anything loads, for run and guard (read has its own): rows K30 and R53, a clean tracked tree, HEAD the
    commit given, the clock before the cutoff (exit 7 at or after it), the frozen T3 files at their pins. Returns
    (row, stamp)."""
    row = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    head_ok(st, commit)
    if past_cutoff():
        print(f"refused: {LATE_REFUSAL}", file=sys.stderr, flush=True)
        sys.exit(EXIT_LATE)
    k22.frozen_ok()
    return row, st


def reading_exists() -> bool:
    """A reading of row K30 is in its folder (RESULTS, looked up at call time)."""
    return any((Path(RESULTS) / n).exists() for n in ("reading.json", "reading.txt"))


def step_passed(out: Path, name: str, st: dict, what: str) -> None:
    """A scoring part starts only after a passing guard.json of this commit and machine (exit 2 otherwise)."""
    p = Path(out) / name
    if not p.is_file():
        refuse(f"{p} is missing: {what} runs and passes first")
    try:
        g = json.loads(p.read_text())
    except ValueError:
        refuse(f"{p} is not JSON")
    if not isinstance(g, dict) or g.get("passes") is not True:
        refuse(f"{what} did not pass: no part runs after it")
    prov = g.get("provenance") or {}
    if prov.get("commit") != st["commit"] or prov.get("machine") != st["machine"]:
        refuse(f"{name} is from another commit or machine: refused")


# ---------------------------------------------------------------------------------------------------------------------
# The depth table and the arms' ratio vectors


def weighted_geomean(ratios: dict[str, float], weights: dict[str, float]) -> float:
    """exp(sum_t w_t log r_t / sum_t w_t) over the labels of `ratios` (every one must have a weight)."""
    w = [float(weights[t]) for t in ratios]
    tot = math.fsum(w)
    if not tot > 0:
        raise Refused("the weights of the geometric mean sum to zero")
    return math.exp(math.fsum(wi * math.log(float(r)) for wi, r in zip(w, ratios.values(), strict=True)) / tot)


def load_table(path: Path | None = None) -> dict:
    """The committed table at its pin (TABLE_SHA256, looked up at call time): per task the other replicate's ratio per
    label, WT E8.75's cells per label and the global ratio over all labels (1.0 where the frozen run measured none).
    Refused when the file is missing or off its pin, a column is missing, a label repeats, a ratio is not finite and
    positive, or a weight is not finite and non-negative or the weights sum to zero."""
    p = TABLE if path is None else Path(path)
    try:
        raw = p.read_bytes()
    except OSError as e:
        raise Refused(f"the depth table {p} cannot be read ({e.__class__.__name__})") from None
    sha = hashlib.sha256(raw).hexdigest()
    if sha != TABLE_SHA256:
        raise Refused(f"the depth table {p.name} is not at its pinned sha256 ({TABLE_SHA256[:12]}...)")
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    need = ["label", WEIGHTS, *COLUMN.values()]
    if not rows or any(c not in rows[0] for c in need):
        raise Refused(f"the depth table lacks a column of {need}")
    labels = [r["label"] for r in rows]
    if len(set(labels)) != len(labels) or any(not t for t in labels):
        raise Refused("the depth table's labels are empty or repeat")

    def num(r: dict, c: str) -> float:
        try:
            return float(r[c])
        except (TypeError, ValueError):
            raise Refused(f"the depth table's {c} for {r['label']!r} is not a number") from None

    weights = {r["label"]: num(r, WEIGHTS) for r in rows}
    if any(not math.isfinite(w) or w < 0 for w in weights.values()) or not math.fsum(weights.values()) > 0:
        raise Refused(f"the depth table's {WEIGHTS} are not finite and non-negative with a positive sum")
    ratios = {}
    for n, c in COLUMN.items():
        ratios[n] = {r["label"]: num(r, c) for r in rows}
        if any(not math.isfinite(x) or x <= 0 for x in ratios[n].values()):
            raise Refused(f"the depth table's {c} holds a ratio that is not finite and positive")
    return {"sha256": sha, "labels": labels, "weights": weights, "ratios": ratios,
            "global": {n: weighted_geomean(ratios[n], weights) for n in TASKS}}


def listed(table: dict, task: str) -> dict[str, float]:
    """The other replicate's ratios that differ from 1.0, as the frozen submission_parameters lists the pooled ones (a
    label at 1.0 is one the frozen run measured no ratio for, or one at no change: either way 1.0 on the vocabulary)."""
    return {t: r for t, r in table["ratios"][task].items() if r != 1.0}


def depth_vectors(bs, V, table: dict) -> dict:
    """Per task and depth arm, the ratio vector on the vocabulary V (one entry per label, in V's order): dx1 the
    builder's depth_ratios of the listed ratios, dxh that ** 0.5, dxg the global ratio everywhere. A listed label that V
    lacks is refused by depth_ratios (build_submission.Refused)."""
    out = {}
    for n in TASKS:
        r1 = np.asarray(bs.depth_ratios(V, {"depth_ratio_pooled": listed(table, n)}), dtype=np.float64)
        vec = {a: r1 ** e for a, e in EXPONENT.items()}
        vec[GLOBAL_ARM] = np.full(V.L, table["global"][n], dtype=np.float64)
        if any(v.shape != (V.L,) or not np.all(np.isfinite(v) & (v > 0)) for v in vec.values()):
            raise RuntimeError("a depth vector is not finite and positive on every label of the vocabulary")
        out[n] = vec
    return out


def depth_record(V, table: dict, vectors: dict, task: str | None = None) -> dict:
    """Target-free: the table's pin, the pairing, the global ratios and the vectors on V by label (one task's, or both
    tasks' when task is None), with the table's labels this vocabulary lacks and whether V's cell counts are the
    table's weights (reported only)."""
    tasks = list(TASKS) if task is None else [task]
    count = dict(zip(V.labels, (int(x) for x in np.asarray(V.count)), strict=True)) if hasattr(V, "count") else {}
    return {"table": TABLE_REL, "table_sha256": table["sha256"], "pairing": dict(COLUMN),
            "global": {n: float(table["global"][n]) for n in TASKS},
            "listed_labels": {n: len(listed(table, n)) for n in TASKS},
            "table_labels_not_in_vocabulary": [t for t in table["labels"] if t not in set(V.labels)],
            "vocabulary_counts_are_the_weights": bool(count) and all(
                count.get(t) == table["weights"][t] for t in table["labels"] if t in count),
            "vectors": {n: {a: {t: float(x) for t, x in zip(V.labels, vectors[n][a], strict=True)}
                            for a in DEPTH_ARMS} for n in tasks}}


def load_depth(path: Path | None = None) -> dict:
    """The table, checked and held in _K30 for fold_factors."""
    table = load_table(path)
    _K30["table"] = table
    return table


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K18's fold_factors (the original, called unchanged) plus the three arms' ratio vectors for both tasks on this
    world's reference vocabulary, and their target-free record."""
    table = _K30.get("table")
    if not table:
        raise RuntimeError("the depth table is not loaded (k30_depth_read.py loads it before any arm)")
    F = _K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    vectors = depth_vectors(bs, ref.V, table)
    F["k30"] = {"ratio": vectors, "V": ref.V, "record": depth_record(ref.V, table, vectors)}
    return F


@contextmanager
def task_set(task: str):
    """The task whose ratios the depth arms take, for the builds inside the context; cleared on exit and on an
    exception, so no build outside it can take a task's ratios."""
    if task not in TASKS:
        raise RuntimeError(f"unknown task {task!r}")
    _K30["task"] = task
    try:
        yield
    finally:
        _K30.pop("task", None)


def _counts_figure(C) -> dict:
    """Mean counts per cell, the number of cells and a sha256 of the counts (their dtype, shape and bytes, hashed in
    place: no copy of a large matrix)."""
    C = np.ascontiguousarray(np.asarray(C))
    h = hashlib.sha256(f"{C.dtype.str}{C.shape}".encode())
    if C.size:
        h.update(memoryview(C).cast("B"))
    return {"mean_counts": float(C.sum(1, dtype=np.int64).mean()) if C.shape[0] else 0.0, "cells": int(C.shape[0]),
            "sha256": h.hexdigest()}


@contextmanager
def observe_counts(rec: dict):
    """For the build inside: the counts t3_decompose.build hands rescale_depth (before the depth step) and lognorm (the
    cells as built), recorded as their mean per cell, their number and their sha256. The frozen functions are called
    unchanged and restored on exit, as build_submission.observe_draws observes draw_rows."""
    t3d = harness._t3_head().t3d
    real_rescale, real_lognorm = t3d.rescale_depth, t3d.lognorm

    def rescale(C, ratio, P, rng):
        rec["pre_step"] = _counts_figure(C)
        return real_rescale(C, ratio, P, rng)

    def lognorm(C):
        rec["built"] = _counts_figure(C)
        return real_lognorm(C)

    t3d.rescale_depth, t3d.lognorm = rescale, lognorm
    try:
        yield rec
    finally:
        t3d.rescale_depth, t3d.lognorm = real_rescale, real_lognorm


def _note_counts(arm: str, seed: int, rec: dict) -> None:
    """Keep one build's figures for the task in progress (none outside task_set), until the seed's task file takes
    them (run_plain and _dm_task start a fresh record for each seed)."""
    if "task" not in _K30 or "built" not in rec:
        return
    fig = {"seed": int(seed), "mean_counts": rec["built"]["mean_counts"], "cells": rec["built"]["cells"],
           "built_sha256": rec["built"]["sha256"]}
    if "pre_step" in rec:
        fig |= {"pre_step_mean_counts": rec["pre_step"]["mean_counts"], "pre_step_sha256": rec["pre_step"]["sha256"]}
    _K30.setdefault("counts", {})[arm] = fig


def counts_record() -> dict:
    """The figures of the builds since task_set began: per arm its mean counts per built cell (and, for a depth arm,
    before the step), and per depth arm whether its counts before the step are byte for byte the base's built
    counts at the same seed (K18's pairing: rows and factor draws shared; reported, deciding nothing)."""
    figs = dict(_K30.get("counts") or {})
    base = figs.get(BASE) or {}
    same = {a: bool(base) and figs[a].get("seed") == base.get("seed")
            and figs[a].get("pre_step_sha256") == base.get("built_sha256") for a in DEPTH_ARMS if a in figs}
    return {"arms": figs, "pre_step_equals_base": same}


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """One arm's drawn cells, its counts observed. The base: K18's own build_arm for barred_form. A depth arm: the same
    call (stable_rng(seed, STREAM), prediction_inputs' keywords from the barred spec, which carry no ratio, no stage
    mixing and exponent 1) with construct's `ratio` set to the arm's vector for the task in progress."""
    if arm not in ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    rec: dict = {}
    if arm == BASE:
        with observe_counts(rec):
            built = _K18_ORIGINAL["build_arm"](bs, ref, k18.BARRED, n, seed, F)
        _note_counts(arm, seed, rec)
        return built
    task = _K30.get("task")
    if task not in TASKS:
        raise RuntimeError("no task is set for the depth step (build through this script's wrappers)")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["barred_spec"], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    k30 = F["k30"]
    if k30["V"] is not ref.V:
        raise RuntimeError("the depth vectors were made on another vocabulary than this reference's")
    kw["ratio"] = k30["ratio"][task][arm]
    with observe_counts(rec):
        built = bs.construct(ref, n, rng, **kw)
    _note_counts(arm, seed, rec)
    return built


_K18_ORIGINAL = {"ARMS": k18.ARMS, "build_arm": k18.build_arm, "fold_factors": k18.fold_factors}


@contextmanager
def k18_arms(arms: tuple[str, ...]):
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's, restored on exit and on an
    exception. K18's file is not changed; its plain build and depth-matched row then build this try's arms."""
    if not set(arms) <= set(ARMS) or BASE not in arms:
        raise RuntimeError(f"arms {arms} are not a set of this try's arms with the base")
    k18.ARMS, k18.build_arm, k18.fold_factors = tuple(arms), build_arm, fold_factors
    k18._FOLD.clear()
    try:
        yield
    finally:
        k18.ARMS = _K18_ORIGINAL["ARMS"]
        k18.build_arm = _K18_ORIGINAL["build_arm"]
        k18.fold_factors = _K18_ORIGINAL["fold_factors"]
        k18._FOLD.clear()


def part_arms(part: str) -> tuple[str, ...]:
    return (BASE,) if part == "guard" else ARMS


TASK_OF_KEY = {t.key: n for n, t in TASKS.items()}
TASK_OF_TARGET = {v: k for k, v in k18.DM_TARGETS.items()}  # rep1 -> g1, rep2 -> g2


def _build_plain(ctx, arm: str):
    """K18's plain build with the task taken from the harness's task key (its count figures join the seed's record,
    which run_plain starts and takes)."""
    with task_set(TASK_OF_KEY[ctx.task.key]):
        return k18._build_plain(ctx, arm)


def _dm_task(task: tuple[int, int]) -> dict:
    """K18's depth-matched row for (target index, seed), the task taken from the target's name, with the row's count
    figures (a fresh record for the row)."""
    ti, _ = task
    _K30["counts"] = {}
    with task_set(TASK_OF_TARGET[k18._STATE["ctx"].targets[ti].name]):
        row = k18._dm_task(task)
    return row | {"counts": counts_record()}


def candidates(arms: tuple[str, ...]) -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build_plain(ctx, arm),
                                 params={"form": arm, "stream": STREAM,
                                         "factor": "the barred s1's (build_submission.py --branch barred --submission "
                                                   "s1)",
                                         "depth": None if arm == BASE else (
                                             "global" if arm == GLOBAL_ARM else f"per label, exponent {EXPONENT[arm]}"),
                                         "pairing": None if arm == BASE else dict(COLUMN), "guard": arm == BASE},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


# ---------------------------------------------------------------------------------------------------------------------
# Running


def progress_log(part: str, world: str, task: str):
    """The harness's per-seed line lists every arm's score; this prints the part (guard only), world, task and seed."""
    lead = "K30 guard" if part == "guard" else "K30"

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"{lead} {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


def k18_rows(refs: Path = K18_TASKS) -> tuple[dict, list[str]]:
    """K18's committed barred_form rows per (world, task), from its task files at this script's pins (K21's reader)."""
    saved = dict(k21.K18_SHA256)
    k21.K18_SHA256.clear()
    k21.K18_SHA256.update(K18_SHA256)
    try:
        return k21.k18_rows(refs)
    finally:
        k21.K18_SHA256.clear()
        k21.K18_SHA256.update(saved)


def guard_refs(refs: Path = K18_TASKS) -> dict:
    """What the guard compares with: K18's committed rows, refused when a file is missing or off its pin."""
    rows, why = k18_rows(refs)
    if why:
        raise Refused("; ".join(why) + ": refused")
    return rows


def mark_run_started(out: Path, st: dict, weight: str) -> None:
    """The run step's start, which counts the three tries: R/run_started.json, made once by the first scoring part
    that passes its refusals (an exclusive create; the other parts find it there)."""
    try:
        with open(Path(out) / RUN_MARK, "x") as fh:
            json.dump({"row": ROW_LABEL, "ruling": RULING, "counted": list(COUNTED), "s3_weight": weight,
                       "commit": st["commit"], "machine": st["machine"],
                       "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}, fh)
    except FileExistsError:
        pass


def run_plain(a, st: dict, out: Path, part: str) -> None:
    """One plain task, seed by seed through harness.run_task (K8's harness), written as
    <root>/plain_<task>/seed_<s>.json with the task's depth record, the seed's count figures and the S3 weight; seeds
    already written at this commit on this machine are skipped."""
    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / f"plain_{a.task}"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    arms = part_arms(part)
    lead = "K30 guard" if part == "guard" else "K30"
    with k18_arms(arms):
        for s in k21.parse_seeds(a.seeds):
            path = folder_ / f"seed_{s}.json"
            if path.exists():
                print(f"{lead} plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            _K30["counts"] = {}
            res = harness.run_task(TASKS[a.task], candidates(arms), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log(part, "plain", a.task))
            got = res["provenance"]
            if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
                raise k22._Unclean(f"the harness's stamp says commit {got.get('commit')}, clean "
                                   f"{got.get('tree_clean')}: no task file written (INCOMPLETE)")
            k22._stamp_ok(st)
            fold = k18._FOLD.get("F")
            if not fold or "k30" not in fold:
                raise RuntimeError("the fold's depth vectors were not made: no task file written")
            res["provenance"]["machine"] = st["machine"]
            res |= {"world": "plain", "part": part, "task_source": SCRIPT, "seed": s, "arms": list(arms),
                    "depth": depth_record(fold["k30"]["V"], _K30["table"], fold["k30"]["ratio"], a.task),
                    "counts": counts_record(), "s3_weight": a.s3_weight, "barred_factor": rec}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)


def u3d_distance(ctx, vectors: dict) -> dict:
    """Per task, the largest |dx1 vector - the live context's other-replicate ratio| (U3d's S1.ratio for that target,
    ctx.targets[other].full.ratio): the table's column and U3d's ratio are one quantity (reported)."""
    out = {}
    for n in TASKS:
        other = next(T for T in ctx.targets if T.name == k18.DM_TARGETS[OTHER_TASK[n]])
        out[n] = float(np.max(np.abs(vectors[n]["dx1"] - np.asarray(other.full.ratio, dtype=np.float64))))
    return out


def run_dm(a, st: dict, out: Path, part: str) -> None:
    """The depth-matched world, set up as k18.run_dm and rows K21 to K24 and K29 set it up (the same functions in the
    same order), then each (replicate, seed) row through k18.dm_row with this part's arms, the task set from the
    target's name; rows already written at this commit on this machine are skipped."""
    import multiprocessing as mp

    import orth_diag as od

    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / "dm"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    seeds = k21.parse_seeds(a.seeds)
    bs = harness._t3_head()
    t3d = bs.t3d
    spec, dec, rules = k18.barred_spec(Path(a.decision))
    if list(rules["protocol"]["seeds"]) != list(SEEDS):
        raise RuntimeError("the frozen script's seeds are not vec.adoption.SEEDS")
    paths = od.context_paths(Path(a.data), with_targets=True)
    n_par = a.workers * a.threads
    arms = part_arms(part)
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=Path(tmp))
        if sorted(T.name for T in ctx.targets) != sorted(k18.DM_TARGETS.values()):
            raise RuntimeError(f"targets {[T.name for T in ctx.targets]}, expected {sorted(k18.DM_TARGETS.values())}")
        k18.quiet_selfcheck(t3d, ctx)
        genes = list(ctx.genes)
        phi_b, frec = k18.barred_factor(Path(a.run_dir), dec, spec, genes)
        ref_C, ref_lab, ref_P, gap = od.reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        wt_path = Path(paths["wt"])
        series = bs.wt_revert_series(wt_path, {n: Path(a.data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                               if n != wt_path.stem})
        phi, per, group, rrec = bs.wt_revert_factors(wt_path, series, genes, ctx.wt.pool.C)
        refw = bs.Reference(preset=ctx.preset_name, rules=rules, genes=genes, V=ctx.V, pool=ctx.wt.pool,
                            coords=np.zeros((ctx.wt.pool.C.shape[0], 3)))
        files = [*(Path(a.data) / f for f in od.WT_FILES.values()), *(Path(a.data) / f for f in od.KO_FILES.values()),
                 *series.values(), Path(a.decision), Path(a.run_dir) / "factors" / f"{od.BARRED}.csv"]
        with k18_arms(arms):
            F = k18.fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b,
                                 spec)
            u3d = u3d_distance(ctx, F["k30"]["ratio"])
            setup = {"world": "dm", "part": part, "barred_factor": frec,
                     "wt_revert": {k: v for k, v in rrec.items() if k != "series"},
                     "reference_counts_largest_gap": gap, "depth": F["k30"]["record"],
                     "u3d_other_ratio_largest_distance": u3d, "s3_weight": a.s3_weight,
                     "files_sha256": {str(f): t3d.sha256(f) for f in dict.fromkeys(files)}}
            sp = folder_ / "setup.json"
            if sp.exists():
                old = json.loads(sp.read_text())
                if {k: old.get(k) for k in setup} != json.loads(json.dumps(setup, default=float)):
                    sys.exit(f"{sp} records another set-up (inputs, factor, depth vectors or weight): not run")
            else:
                k18.write_task_file(sp, setup | {"provenance": st})
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder_ / f"{TASK_OF_TARGET[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                k22._stamp_ok(st)
                task = TASK_OF_TARGET[row["target"]]
                path = folder_ / f"{task}_seed_{row['seed']}.json"
                k18.write_task_file(path, row | {"world": "dm", "part": part, "task": task, "arms": list(arms),
                                                 "depth": depth_record(ctx.V, _K30["table"], F["k30"]["ratio"], task),
                                                 "u3d_other_ratio_largest_distance": u3d[task],
                                                 "s3_weight": a.s3_weight, "provenance": st, "task_source": SCRIPT})
                print(f"wrote {path}", flush=True)

            if a.workers > 1 and len(todo) > 1:
                if not bs.openmp_limited():
                    raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
                with mp.get_context("fork").Pool(min(a.workers, len(todo))) as p:
                    for row in p.imap_unordered(_dm_task, todo):
                        write(row)
            else:
                for t in todo:
                    write(_dm_task(t))


def cmd_run(a) -> None:
    if a.world == "plain" and a.task not in TASKS:
        refuse("--world plain needs --task g1 or g2")
    if a.world != "plain" and a.task is not None:
        refuse("--task is for --world plain only")
    row, st = preconditions(a.commit)
    weight_ok(row, a.s3_weight)
    if reading_exists():  # as read: once a reading is there, no arm is scored again, in this folder or another
        print(f"refused: {RESULTS} holds a reading already; no arm is scored again", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    out = Path(a.out)
    try:
        if a.part == "guard":
            guard_refs()
        else:
            step_passed(out, GUARD_FILE, st, "the guard part")
        load_depth()
        out.mkdir(parents=True, exist_ok=True)
        if a.part == "score":
            mark_run_started(out, st, a.s3_weight)
        (run_plain if a.world == "plain" else run_dm)(a, st, out, a.part)
    except Refused as e:
        refuse(str(e))
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seeds written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)


# ---------------------------------------------------------------------------------------------------------------------
# The guard


def guard_files() -> list[str]:
    return ([f"{GUARD_DIR}/plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"{GUARD_DIR}/dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def score_files() -> list[str]:
    return ([f"plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def cell_of(rel: str) -> tuple[str, int]:
    """(cell, seed) of a task file's path relative to its part's root."""
    seed = int(rel.rsplit("_", 1)[1][:-len(".json")])
    if rel.startswith("plain_"):
        return f"plain {rel.split('/')[0][len('plain_'):]}", seed
    if rel.startswith("dm/"):
        return f"dm {rel[3:5]}", seed
    raise ValueError(f"unknown task file {rel}")


def rows_of(runs: Path, files: list[str], prefix: str, part: str, want: tuple[str, ...], why: list[str], sha: dict,
            prov: dict, depth: dict | None = None, weights: dict | None = None) -> dict:
    """Rows per cell from task files (plain: rows[0] with the file's count figures; dm: the row), each checked for its
    part, task, seed and arms (the plain ones the floor too); a file that fails is a reason in `why` and gives no row.
    Each file's depth record and S3 weight go into `depth` and `weights` when given."""
    rows: dict = {c: [] for c in CELLS}
    for f in files:
        r, sha[f] = k21._read_json(Path(runs) / f)
        prov[f] = r.get("provenance") or {}
        if depth is not None:
            depth[f] = r.get("depth")
        if weights is not None:
            weights[f] = r.get("s3_weight")
        cell, seed = cell_of(f[len(prefix):])
        world, task = cell.split()
        if r.get("part") != part:
            why.append(f"{f}: not a file of the {part} part")
            continue
        if world == "plain":
            if (r.get("task") or {}).get("key") != TASKS[task].key or [q.get("seed") for q in r.get("rows") or []] \
                    != [seed]:
                why.append(f"{f}: another task or seed")
                continue
            row, need = {**r["rows"][0], "counts": r.get("counts")}, (harness.FLOOR, *want)
        else:
            if r.get("task") != task or r.get("seed") != seed:
                why.append(f"{f}: another task or seed")
                continue
            row, need = r, want
        miss = [x for x in need if x not in (row.get("scores") or {}) or x not in (row.get("skills") or {})
                or x not in (row.get("t3") or {})]
        if miss:
            why.append(f"{f}: missing arms {miss}")
            continue
        rows[cell].append(row)
    return {c: sorted(v, key=lambda q: q["seed"]) for c, v in rows.items()}


def guard_eval(rows: dict, refs: dict) -> dict:
    """The guard: the base within K17's tolerance (k18.guard_passes: ten-seed mean within 0.1, every seed within 0.5)
    of K18's committed barred_form in all four cells (the ten seeds needed in each)."""
    cells = {}
    for c in CELLS:
        rr = rows.get(c) or []
        ok = [q["seed"] for q in rr] == list(SEEDS)
        d = k18.distance(rr, BASE, refs.get(tuple(c.split())) or [], k18.BARRED) if ok else None
        cells[c] = {"distance": d, "passes": bool(ok and k18.guard_passes(d))}
    passes = all(v["passes"] for v in cells.values())
    return {"cells": cells, "failed": [] if passes else [BASE], "passes": passes}


def match_eval(rows: dict, grows: dict) -> dict:
    """The run step's base against the guard step's, seed by seed, within the same tolerance in all four cells."""
    cells = {}
    for c in CELLS:
        rr, gg = rows.get(c) or [], grows.get(c) or []
        ok = [q["seed"] for q in rr] == list(SEEDS) and [q["seed"] for q in gg] == list(SEEDS)
        d = k18.distance(rr, BASE, gg, BASE) if ok else None
        cells[c] = {"distance": d, "passes": bool(ok and k18.guard_passes(d))}
    return {"cells": cells, "passes": all(v["passes"] for v in cells.values())}


def cmd_guard(a) -> None:
    """The guard check over R/guard/: R/guard.json (made once); exit 0 when the guard holds, 6 when it fails (VOID: no
    scoring part starts; the base counts one try), 3 when a guard file is missing, 2 when the files are not this
    commit's and machine's or K18's files are off their pins, 7 at or after the cutoff."""
    row, st = preconditions(a.commit)
    runs = Path(a.runs)
    p = runs / GUARD_FILE
    if p.exists():
        print(f"refused: {p} exists: the guard check is made once", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    missing = [f for f in guard_files() if not (runs / f).exists()]
    if missing:
        print(f"INCOMPLETE: {len(missing)} guard file(s) missing", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    try:
        refs = guard_refs()
    except Refused as e:
        refuse(str(e))
    try:
        why: list[str] = []
        sha, prov = {}, {}
        rows = rows_of(runs, guard_files(), f"{GUARD_DIR}/", "guard", (BASE,), why, sha, prov)
        bad = [f for f, q in prov.items() if q.get("tree_clean") is not True or q.get("commit") != st["commit"]
               or q.get("machine") != st["machine"]]
        if bad or why:
            refuse(f"the guard files are not this commit's and machine's, or not well formed ({(bad + why)[:4]})")
        res = guard_eval(rows, refs)
        k18.write_task_file(p, {**res, "row": row["id"], "files_sha256": sha, "provenance": st})
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print(f"guard: {'passes' if res['passes'] else 'FAILS'}", flush=True)
    if not res["passes"]:
        sys.exit(EXIT_GUARD)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def counts_summary(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """Per arm (not the floor): its mean counts per built cell over the seeds, against barred_form's, and the seeds
    whose counts before the step are barred_form's built counts (reported; None where a seed has no figure)."""
    def mean_counts(q: dict, a: str):
        return (((q.get("counts") or {}).get("arms") or {}).get(a) or {}).get("mean_counts")

    base = [mean_counts(q, BASE) for q in rows]
    mb = float(np.mean(base)) if base and all(x is not None for x in base) else None
    out = {}
    for a in arms:
        if a == harness.FLOOR:
            continue
        v = [mean_counts(q, a) for q in rows]
        if not v or any(x is None for x in v):
            out[a] = None
            continue
        m = float(np.mean(v))
        same = None if a == BASE else sum(bool(((q.get("counts") or {}).get("pre_step_equals_base") or {}).get(a))
                                          for q in rows)
        out[a] = {"mean_counts": m, "against_base": m / mb if mb else None, "pre_step_equals_base_seeds": same}
    return out


def summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """Per arm: the score's mean and sd, each metric's mean skill, the severity states, seeds scored and gated, the
    median beta and r2. Paired (k7._paired): each depth arm against the base (D_a), dx1 - dxh, dx1 - dxg, the base
    against the floor (plain). Per-term: each pair's mean per-seed skill difference. Gate flips: per depth arm, the
    seeds where its severity state is not the base's. Counts: counts_summary."""
    out: dict = {"arms": {}, "paired": {}, "per_term": {}, "flips": {}, "counts": counts_summary(rows, arms)}
    for x in arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "severity_states": {v: sum(z["state"] == v for z in sev) for v in sorted({z["state"] for z in sev})},
            "severity_scored": int(sum(v["state"] == "scored" for v in sev)), "severity_gated": k21._gated(sev),
            "median_beta": k18._median([v.get("beta") for v in sev]),
            "median_r2": k18._median([v.get("r2") for v in sev])}
    pairs = {f"{a} - {BASE}": (a, BASE) for a in DEPTH_ARMS} | {f"{x} - {y}": (x, y) for x, y in REPORTED_PAIRS}
    if harness.FLOOR in arms:
        pairs[f"{BASE} - {harness.FLOOR}"] = (BASE, harness.FLOOR)
    out["paired"] = {k: k7._paired(rows, x, y) for k, (x, y) in pairs.items()}
    for k, (x, y) in pairs.items():
        if y == harness.FLOOR:
            continue
        out["per_term"][k] = {m: float(np.mean([q["skills"][x][m] - q["skills"][y][m] for q in rows]))
                              for m in rows[0]["skills"][y]}
    for a in DEPTH_ARMS:
        out["flips"][a] = [q["seed"] for q in rows
                           if q["t3"][a]["severity"]["state"] != q["t3"][BASE]["severity"]["state"]]
    return out


def ahead(v: dict) -> bool:
    """AHEAD: a paired mean above +1.0 with at least 8 of 10 seeds positive."""
    return bool(v["mean"] > NAME_GAIN and v["ahead"] >= MIN_AHEAD)


def _avg(xs: list[float]) -> Fraction:
    return sum((Fraction(float(x)) for x in xs), Fraction(0)) / len(xs)


def ev(plain_means: list[float], dm_means: list[float], weight: str) -> Fraction:
    """EV = w x mean(plain cells) + (1 - w) x mean(dm cells), w the S3 weight (1/2 for 0.5, 1 for plain), exact in
    fractions of the float means."""
    w = S3_WEIGHTS[weight]
    return w * _avg(plain_means) + (1 - w) * _avg(dm_means)


def decide(per: dict, weight: str) -> dict:
    """The rule fixed in row K30 on summarise()'s output per cell: AHEAD, EV at the S3 weight, S1 to S4 and NAME
    (module docstring)."""
    if weight not in S3_WEIGHTS:
        raise ValueError(f"unknown S3 weight {weight!r}")
    d = {a: {c: per[c]["paired"][f"{a} - {BASE}"] for c in CELLS} for a in DEPTH_ARMS}
    up = {a: {c: ahead(v) for c, v in d[a].items()} for a in DEPTH_ARMS}
    plain_mean = {a: [d[a][c]["mean"] for c in PLAIN_CELLS] for a in DEPTH_ARMS}
    dm_mean = {a: [d[a][c]["mean"] for c in DM_CELLS] for a in DEPTH_ARMS}
    evs = {a: ev(plain_mean[a], dm_mean[a], weight) for a in DEPTH_ARMS}
    gated = {c: {x: int(per[c]["arms"][x]["severity_gated"]) for x in ARMS} for c in CELLS}
    pass_s1 = {a: all(up[a][c] for c in PLAIN_CELLS) for a in DEPTH_ARMS}
    over = {a: [c for c in PLAIN_CELLS if gated[c][a] > gated[c][BASE] + MAX_EXTRA_GATED] for a in DEPTH_ARMS}
    pass_s2 = {a: not over[a] for a in DEPTH_ARMS}
    may = [a for a in DEPTH_ARMS if pass_s1[a] and pass_s2[a]]
    order = {a: i for i, a in enumerate(DEPTH_ARMS)}
    best = max(may, key=lambda a: (evs[a], -order[a])) if may else None  # an exact tie: the first of dx1, dxh, dxg
    stop = {"S1": not any(pass_s1.values()),
            "S2": any(pass_s1.values()) and not may,
            "S3": best is not None and evs[best] <= EV_BAR,
            "S4": all(abs(d[a][c]["mean"]) < NULL_BAND for a in DEPTH_ARMS for c in PLAIN_CELLS)}
    clauses = [k for k, v in stop.items() if v]
    named = None if clauses else best
    ties = sorted(a for a in may if a != best and evs[a] == evs[best]) if best else []
    return {"stop": stop, "clauses": clauses, "named": named, "s3_weight": weight,
            "outcome": STOP_LINE.format(", ".join(clauses)) if clauses else NAME_LINE.format(arm=named),
            "file": FILE_OF[named] if named else None,
            "arms": {a: {"ev": float(evs[a]), "ev_exact": f"{evs[a].numerator}/{evs[a].denominator}",
                         "plain_mean": float(_avg(plain_mean[a])), "dm_mean": float(_avg(dm_mean[a])),
                         "passes_s1": pass_s1[a], "passes_s2": pass_s2[a], "s2_cells": over[a],
                         "may_be_named": a in may} for a in DEPTH_ARMS},
            "best": best, "ev_tied_with_best": ties,
            "table": {a: {c: {"mean": float(v["mean"]), "sd": float(v["sd"]), "ahead_seeds": int(v["ahead"]),
                              "n": int(v["n"]), "AHEAD": up[a][c]} for c, v in d[a].items()} for a in DEPTH_ARMS},
            "gated": gated}


INCOMPLETE = ("Try K30: INCOMPLETE: not a reading; guard.json, a guard file or a seed file is missing (a crash, a "
              "refusal such as a dirty tree, out of memory or a kill); nothing was written. If the cause is outside "
              "the code, rerun the missing part into this folder at the run's commit on the same machine (the runs "
              "skip what is there); after a committed code fix, rerun every step at the new commit into a fresh "
              f"folder and never read this one. At or after {CUTOFF_TEXT} a read writes the VOID reading '{LATE}' "
              "instead")


def run_started(runs: Path) -> bool:
    """The run step has started: R/run_started.json, or a scoring file, is there."""
    return (Path(runs) / RUN_MARK).exists() or any((Path(runs) / f).exists() for f in score_files())


def _missing(missing: list[str], counted: list[str], started: bool, base: dict) -> dict:
    """Files missing: INCOMPLETE (not a reading) before the cutoff; at or after it, the VOID reading 'INCOMPLETE past
    its cutoff', with the tries already counted."""
    if not past_cutoff():
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    why = [f"{LATE} ({CUTOFF_TEXT}): {len(missing)} file(s) missing, e.g. {', '.join(missing[:4])}"]
    return base | {"void": True, "late": True, "missing": missing, "why": why, "tries_counted": len(counted),
                   "counted_arms": counted, "run_step_started": started, "outcome": VOID_LINE.format(why[0])}


def expected_depth(table: dict) -> dict:
    """What every task file's depth record must say, vocabulary aside: the pin, the pairing, the global ratios, and per
    task and arm each label's ratio (dx1 the listed ratio or 1.0, dxh its square root, dxg the global ratio)."""
    def want(task: str, arm: str, label: str) -> float:
        if arm == GLOBAL_ARM:
            return float(table["global"][task])
        return float(np.float64(listed(table, task).get(label, 1.0)) ** EXPONENT[arm])

    return {"sha256": table["sha256"], "pairing": dict(COLUMN), "global": {n: float(table["global"][n]) for n in TASKS},
            "want": want}


def depth_why(depth: dict, files: list[str], table: dict | None) -> list[str]:
    """Reasons a depth record is not the pinned table's (empty: every record is): its pin, pairing and global ratios
    exactly, and every label's ratio of its own task's three vectors within a relative 1e-12."""
    if table is None:
        return ["the depth table is not at its pinned sha256 or not well formed at read time"]
    exp = expected_depth(table)
    why = []
    for f in files:
        r = depth.get(f)
        if not isinstance(r, dict) or r.get("table_sha256") != exp["sha256"] or r.get("pairing") != exp["pairing"] \
                or r.get("global") != exp["global"]:
            why.append(f"{f}: its depth record is not the pinned table's")
            continue
        task = cell_of(f.split("/", 1)[1] if f.startswith(f"{GUARD_DIR}/") else f)[0].split()[1]
        vec = (r.get("vectors") or {}).get(task) or {}
        if set(vec) != set(DEPTH_ARMS) or any(
                not vec[a] or any(not math.isclose(float(x), exp["want"](task, a, t), rel_tol=DEPTH_REL, abs_tol=0.0)
                                  for t, x in vec[a].items()) for a in DEPTH_ARMS):
            why.append(f"{f}: its depth vectors are not the other replicate's ratios of the pinned table")
    return why


def reading(runs: Path, refs: Path = K18_TASKS, *, weight: str) -> dict:
    """Apply the reading of row K30 at the S3 weight to a run folder: the guard step, then the run step. VOID is checked
    first and wins; the STOP and NAME rule is applied only to a run that is not void. No seed file is opened while one
    the reading needs is missing."""
    runs = Path(runs)
    started = run_started(runs)
    sha, prov = {}, {}
    base: dict = {"incomplete": False, "files_sha256": sha, "provenance": prov, "s3_weight": weight}
    three = list(COUNTED) if started else []
    missing = [f for f in (GUARD_FILE, *guard_files()) if not (runs / f).exists()]
    if missing:  # guard.json is opened only once every guard file is there
        return _missing(missing, three, started, base)
    g, sha[GUARD_FILE] = k21._read_json(runs / GUARD_FILE)
    prov[GUARD_FILE] = g.get("provenance") or {}
    base["guard_part"] = {k: g.get(k) for k in ("cells", "failed", "passes")}
    guard_step = g.get("passes") is True
    missing = [f for f in score_files() if not (runs / f).exists()] if guard_step else []
    if missing:
        return _missing(missing, three, started, base)
    why: list[str] = []
    gsha: dict = {}
    depth: dict = {}
    weights: dict = {}
    grows = rows_of(runs, guard_files(), f"{GUARD_DIR}/", "guard", (BASE,), why, gsha, prov, depth, weights)
    if gsha != g.get("files_sha256"):
        why.append("guard.json: its files' sha256 are not the guard files'")
    sha |= gsha
    if not guard_step:  # the run step never starts after a failed guard: one try, the base, whatever else is there
        counted = [BASE]
        why = ["the guard failed in the guard step: barred_form is off K18's committed barred_form", *why]
        return base | {"void": True, "why": why, "tries_counted": len(counted), "counted_arms": counted,
                       "run_step_started": started, "outcome": VOID_LINE.format("; ".join(why))}
    rows = rows_of(runs, score_files(), "", "score", ARMS, why, sha, prov, depth, weights)
    for f, p in prov.items():
        if p.get("tree_clean") is not True:
            why.append(f"{f}: a dirty tree")
    commits = {p.get("commit") for p in prov.values()}
    machines = {p.get("machine") for p in prov.values()}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(map(str, commits))}")
    if len(machines) > 1:
        why.append(f"results from different machines {sorted(map(str, machines))}")
    other = [f for f, w in weights.items() if w != weight]
    if other:
        why.append(f"{len(other)} task file(s) were run under another S3 weight than {weight}, e.g. {other[0]}")
    try:
        table = load_table()
    except Refused:
        table = None
    why += depth_why(depth, [*guard_files(), *score_files()], table)
    k18r, kwhy = k18_rows(refs)
    why += kwhy
    again = None if kwhy else guard_eval(grows, k18r)
    if again is not None and not again["passes"]:
        why.append("the guard: the guard step's barred_form, read again, is off K18's committed barred_form")
    match = match_eval(rows, grows)
    for c, v in match["cells"].items():
        if not v["passes"]:
            dd = v["distance"]
            why.append(f"the run step's barred_form is off the guard step's in {c}" + (
                f" by {dd['mean_diff']:+.3f} in the mean and {dd['max_abs_diff']:.3f} in one seed" if dd else
                " (a seed is missing)") + f" (limits {GUARD_MEAN} and {GUARD_SEED})")
    per: dict = {}
    for c, rr in rows.items():
        if [q["seed"] for q in rr] != list(SEEDS):
            why.append(f"{c}: not the ten seeds")
            continue
        per[c] = summarise(rr, (harness.FLOOR, *ARMS) if c.startswith("plain") else ARMS)
    # the base counts one try when it is measured outside its guard; a cell with no measured distance (a malformed
    # file drops a seed) adds none
    off = [v for m in (again, match) if m for v in m["cells"].values() if v["distance"] is not None and not v["passes"]]
    counted = [*three, *((BASE,) if off else ())]
    mark = runs / RUN_MARK
    u3d = {c: max((float(q["u3d_other_ratio_largest_distance"]) for q in rr
                   if q.get("u3d_other_ratio_largest_distance") is not None), default=None)
           for c, rr in rows.items() if c.startswith("dm")}
    base |= {"guard": again, "guard_match": match, "tries_counted": len(counted), "counted_arms": counted,
             "run_step_started": True, "run_started": json.loads(mark.read_text()) if mark.exists() else None,
             "cells": per, "u3d_other_ratio_largest_distance": u3d,
             "depth": next((r for f, r in depth.items() if f in score_files() and r), None)}
    if why:
        return base | {"void": True, "why": why, "outcome": VOID_LINE.format("; ".join(why))}
    dec = decide(per, weight)
    return base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                   "decision": dec, "outcome": dec["outcome"]}


def reading_code() -> dict:
    st = k18.stamp()
    h = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731
    return {"script": SCRIPT, "script_sha256": h(__file__), "k18_sha256": h(HERE / "k18_built_files.py"),
            "k21_sha256": h(HERE / "k21_factor_forms.py"), "k22_sha256": h(HERE / "k22_outside_forms.py"),
            "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}


def _f(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def _dist(d: dict | None) -> str:
    if d is None:
        return "n/a"
    return (f"mean {d['mean_diff']:+.4f}, largest seed {d['max_abs_diff']:.4f}, seeds differing "
            f"{d['seeds_off'] or 'none'}")


def table_lines(dec: dict) -> list[str]:
    """The decision table: per arm D_a in each cell, EV, S1, S2 and whether it may be named; the severity-gated seeds
    of every arm in every cell."""
    w = dec["s3_weight"]
    lines = [f"per cell, paired by seed, D = arm - barred_form (AHEAD: mean > +1.0 with >= 8 of 10 positive); EV at "
             f"the S3 weight {w} = " + ("0.5 x mean(plain) + 0.5 x mean(dm):" if w == "0.5" else
                                        "mean(plain), the dm cells deciding nothing:")]
    for a in DEPTH_ARMS:
        v = dec["arms"][a]
        lines.append(f"  {a}: " + "; ".join(
            f"{c} {t['mean']:+.2f} (sd {t['sd']:.2f}, {t['ahead_seeds']} of {t['n']} positive)"
            + (" AHEAD" if t["AHEAD"] else "") for c, t in dec["table"][a].items())
            + f"; EV {v['ev']:+.3f} (plain {v['plain_mean']:+.2f}, dm {v['dm_mean']:+.2f}); passes S1 "
            f"{'yes' if v['passes_s1'] else 'no'}, S2 {'yes' if v['passes_s2'] else 'no'}"
            + (f" (more gated than barred_form plus 1 in {', '.join(v['s2_cells'])})" if v["s2_cells"] else "")
            + f"; may be named {'yes' if v['may_be_named'] else 'no'}")
    lines.append("  severity-gated seeds of 10 per arm and cell: " + "; ".join(
        f"{c} " + ", ".join(f"{x} {n}" for x, n in dec["gated"][c].items()) for c in CELLS))
    return lines


def summary_lines(res: dict) -> list[str]:
    """What `read` prints: the outcome line, then for a decided reading the stop clauses and the decision table."""
    lines = [res["outcome"]]
    if not res.get("incomplete") and not res.get("void"):
        dec = res["decision"]
        lines.append("stop: " + ", ".join(f"{k} {'yes' if v else 'no'}" for k, v in dec["stop"].items()))
        lines += table_lines(dec)
    return lines


def text(res: dict) -> str:
    lines = [res["outcome"]]
    if res.get("incomplete"):
        miss = res["missing"]
        return "\n".join([*lines, f"missing: {len(miss)} file(s), e.g. " + ", ".join(miss[:6])]) + "\n"
    lines.append(f"tries counted: {res['tries_counted']} ({', '.join(res['counted_arms']) or 'none'}); the run step "
                 f"started: {res['run_step_started']}; S3 weight {res['s3_weight']}; ruling row {RULING}")
    gp = res.get("guard_part") or {}
    if gp:
        lines.append(f"guard step: passes {gp.get('passes')}, failed {gp.get('failed') or 'none'}")
        for c, v in (gp.get("cells") or {}).items():
            lines.append(f"  guard step {c}: passes {v['passes']}; {_dist(v.get('distance'))}")
    for c, v in ((res.get("guard") or {}).get("cells") or {}).items():
        lines.append(f"  the guard step's barred_form read again, {c}: passes {v['passes']}; {_dist(v['distance'])}")
    for c, v in ((res.get("guard_match") or {}).get("cells") or {}).items():
        lines.append(f"  the run step's barred_form against the guard step's, {c}: passes {v['passes']}; "
                     f"{_dist(v['distance'])}")
    if res["void"]:
        lines.append("VOID (checked first; no STOP or NAME clause is applied): " + "; ".join(res["why"]))
    else:
        dec = res["decision"]
        lines.append(f"commit {res['commit']}, machine {res['machine']}.")
        lines.append("stop: " + ", ".join(f"{k} {'yes' if v else 'no'}" for k, v in dec["stop"].items())
                     + "".join(f"; {CLAUSES[k]}" for k in dec["clauses"]))
        if dec["named"]:
            a = dec["named"]
            v = dec["arms"][a]
            g = dec["gated"]
            lines.append(f"named: {a} (R53 names its file {dec['file']}), EV {v['ev']:+.3f} at the S3 weight "
                         f"{dec['s3_weight']}; the plain cells' mean {v['plain_mean']:+.2f} (the knockout at its run's "
                         f"depth); the dm cells' mean {v['dm_mean']:+.2f} (at its reference's depth) with severity "
                         f"gated in {g['dm g1'][a]} and {g['dm g2'][a]} of 10 seeds (barred_form {g['dm g1'][BASE]} "
                         f"and {g['dm g2'][BASE]}); EV ties with it: {', '.join(dec['ev_tied_with_best']) or 'none'}")
        lines += table_lines(dec)
    lines.append("reported, deciding nothing (each arm's figures, dx1 - dxh, dx1 - dxg, the base against the floor, "
                 "per-term skills, gate flips against the base, mean counts per built cell against barred_form's):")
    for c, v in (res.get("cells") or {}).items():
        lines.append(f"{c}:")
        for a, s in v["arms"].items():
            lines.append(f"    {a:<11} {s['mean']:.2f} (sd {s['sd']:.2f}); severity scored {s['severity_scored']}, "
                         f"gated {s['severity_gated']} {s['severity_states']}, median beta {_f(s['median_beta'])}, "
                         f"r2 {_f(s['median_r2'])}; skills " + ", ".join(f"{m} {x:.3f}" for m, x in
                                                                           s["skills"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
        for k, terms in v["per_term"].items():
            lines.append(f"    per-term skill, {k}: " + ", ".join(f"{m} {x:+.3f}" for m, x in terms.items()))
        lines.append("    gate flips against the base (seeds): " + "; ".join(
            f"{a} {s or 'none'}" for a, s in v["flips"].items()))
        lines.append("    mean counts per built cell (against barred_form's; seeds whose counts before the step are "
                     "barred_form's): " + "; ".join(
                         f"{a} n/a" if s is None else f"{a} {s['mean_counts']:.1f} ({_f(s['against_base'])}"
                         + ("" if s["pre_step_equals_base_seeds"] is None else
                            f"; {s['pre_step_equals_base_seeds']} of {len(SEEDS)}") + ")"
                         for a, s in v["counts"].items()))
    dep = res.get("depth")
    if dep:
        lines.append(f"depth record (target-free): {dep['table']} sha256 {dep['table_sha256'][:12]}; pairing "
                     + ", ".join(f"{n} takes {c}" for n, c in dep["pairing"].items()) + "; global ratios (dxg) "
                     + ", ".join(f"{n} {x:.4f}" for n, x in dep["global"].items()) + "; labels listed "
                     + ", ".join(f"{n} {x}" for n, x in dep["listed_labels"].items()) + "; largest distance of dx1 to "
                     "U3d's other-replicate ratio " + ", ".join(
                         f"{c} {_f(x)}" for c, x in (res.get("u3d_other_ratio_largest_distance") or {}).items()))
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    lines.append("g1 and g2 are two replicates of one knockout run; barred_form is the base and the guard (K18's).")
    if res.get("reading_code"):
        rc = res["reading_code"]
        lines.append(f"read by {rc['script']} {rc['script_sha256'][:12]} (k18 {rc['k18_sha256'][:12]}, k21 "
                     f"{rc['k21_sha256'][:12]}, k22 {rc['k22_sha256'][:12]}) at commit {rc['commit']}, clean "
                     f"{rc['tree_clean']}, machine {rc['machine']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    row = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    head_ok(st, a.commit)
    weight_ok(row, a.s3_weight)
    out = Path(a.out)
    targets = [out / "reading.json", out / "reading.txt"]
    if any(p.exists() for p in targets):
        print(f"refused: {out} holds a reading already; it is never overwritten", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    try:
        res = reading(Path(a.runs), weight=a.s3_weight)
        if not res.get("incomplete"):
            res["reading_code"] = reading_code()
            res["row"] = row["id"]
            res["ruling"] = row["ruling"]
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
        shown = summary_lines(res)
    except Exception:  # noqa: BLE001  (a crash after the files were opened: masked traceback, nothing written)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if res.get("incomplete"):
        print("INCOMPLETE: not a reading; nothing written", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    out.mkdir(parents=True, exist_ok=True)
    for p, body in zip(targets, (js, t), strict=True):
        with open(p, "x") as fh:  # exclusive: an existing reading is never overwritten
            fh.write(body)
    print("\n".join(shown), flush=True)
    print(f"wrote {targets[0]} and {targets[1]}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("guard", "score"), default="score")
    r.add_argument("--world", choices=WORLDS, required=True)
    r.add_argument("--task", choices=sorted(TASKS))
    r.add_argument("--seeds", default=None)
    r.add_argument("--out", required=True)
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--threads", type=int, default=8)
    r.add_argument("--run-dir", default=str(k18.RUN_DIR))
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--decision", default=str(k18.DECISION))
    gd = sub.add_parser("guard")
    gd.add_argument("--runs", required=True)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", default=str(RESULTS))
    for q in (r, d):
        q.add_argument(WEIGHT_OPTION, dest="s3_weight", choices=sorted(S3_WEIGHTS), required=True)
    for q in (r, gd, d):
        q.add_argument("--commit", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "guard": cmd_guard, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
