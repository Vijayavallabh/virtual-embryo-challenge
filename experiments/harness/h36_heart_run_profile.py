#!/usr/bin/env python
"""Row H36 of experiments/harness/TRIES.md: the heart's XPC line (the step-change round's HEART MODEL design, P1, with
every fix of its cold challenge; ruling row R56, the captain's named exception for file OXP125 and his A15 exception,
bound to this read). A labelled diagnostic on one truth. It builds no test file; its last step writes the profile the
builder's `--run-profile <sha256>` takes, and only when everything before it passed.

OXP125 is OTS125 (row R49's file D at trend strength 1.25, the banked heart file) with E10.5's own MERFISH run profile
removed: every drawn cell re-measured under the profile the other runs share, measured against the organisers'
Multiome RNA at the stages both instruments hold (experiments/heart_profile/xpc_table.py) and applied in count space
after the trend step, before the size step (vec.t2_build.apply_run_profile, stream [build seed, 37]).

Rows H28's, H29b's and H34's scripts are imported, not copied: the transition with D's pool as its reference input,
D_f's size and half, the trend step and its scoring are theirs.

Steps, in this order, each refused unless the one before it wrote its verdict at the run commit, none twice (the
launcher, experiments/harness/run_h36.sh, holds "none twice"):
  entry     (E) not a try (row H29b's caveat: it looks at E10.5's labels): every category (cm, endocardium,
            progenitor) at least 200 cells on both instruments at E8.5, E9.5 and E10.5 (E10.5's RNA library 3).
            Exit 0 ENTER, 3 STOP
  profile   (Q) the practice profile, one T1 reading by the strict reading of rows N59 to N61; it reads no E10.5
            file: delta_p = D_E9.5 - D_E8.5, its halves at seed 20261125; STOP unless SD(delta_p) >= 2 sigma_n,
            Spearman of the halves >= 0.50 and the median pairwise Spearman of the delta_p,k >= 0.40. Exit 0 PASS,
            3 STOP
  gate      counting nothing, on `T2h/<=E9.5->E10.5`, seeds 0 to 9: the copy D_f (row H29b's organ:f) and row H34's
            trend_s1.0 and trend_s1.5 on it, rebuilt and scored. G1: D_f reproduces row H29b's committed record on
            this transition exactly ((R) its record, (D) its draw, (M) six plain ten-seed means, (P) 90 paired
            statistics against row H27's committed controls). G2: trend_s1.0's and trend_s1.5's ten-seed means in the
            worlds own and CLEAN equal those of row H34's committed parts (their sha256 pinned here) exactly. Exit 0
            PASS, 3 FAIL (the outcome is then VOID)
  run       the eleven counted arms (X, anti-X, P1 to P9) and their control C, per seed: C is OTS125's practice form,
            h28.trend_expression(ctx, 1.25) on D_f's context, then D_f's size and half; each arm is C's counts with
            one apply_run_profile step (X exp(-delta_p), anti-X exp(+delta_p), P_i delta_p deranged within ten level
            bins of the eligible genes at seeds 20261251 to 20261259), then C's size and half. Eleven tries, counted
            from this step's start
  read      the void conditions first, then the bars, in both worlds (own, CLEAN): B1 X - C >= +0.30 with 9 of 10
            seeds above 0; B2 X - mean(P) >= max(+0.30, 2 x SD of the nine placebos' ten-seed means), 9 of 10 seeds
            above 0, and X's mean above every P_i's; B3 anti-X - mean(P) <= -0.30 with 9 of 10 seeds below 0; B4
            the de_score and de_direction skill differences X - mean(P) both above 0. PASSED AT PRACTICE only if all
            four hold in both worlds; else NOT PASSED
  test-profile  (T) only after the read's PASSED AT PRACTICE is committed (its reading tracked and unchanged at HEAD,
            and the row's outcome cell holding the statement word for word), at that later commit: one heart reading
            (E10.5) and one T1 reading (library 3): delta = D_E10.5 - (D_E8.5 + D_E9.5) / 2 with Q's checks, and
            Spearman against libraries 1 to 3 >= 0.70. PASS writes experiments/heart_profile/results/xpc_profile.json
            and states its sha256; any failure: no file
  refusals  one step's refusals and nothing else (the launcher asks before it marks a step started)
Every step is refused before anything loads unless row H36 is in TRIES.md and well formed (h36_row.py), vec is this
checkout's own, the tracked tree is clean and HEAD is the commit given.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2]  # this checkout
sys.path.insert(0, str(HERE / "src"))

from vec import harness, t2_build  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h34 = _load("h34_heart_organ_trend", HERE / "experiments" / "harness" / "h34_heart_organ_trend.py")
h28, h29b, h26 = h34.h28, h34.h29b, h34.h26
h36_row = _load("h36_row", HERE / "experiments" / "harness" / "h36_row.py")
xpc = _load("xpc_table", HERE / "experiments" / "heart_profile" / "xpc_table.py")

ROW = h36_row.ROW
RULING_ROW, FILE = "R56", h36_row.FILE
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h36-heart-run-profile"
PART = "h36.json"
SOURCE = "experiments/harness/h36_heart_run_profile.py (row H36)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h36_heart_run_profile"
PROFILE = xpc.PROFILE  # what the builder's --run-profile reads
TASK = "T2h/<=E9.5->E10.5"  # the one transition that can host the read (T1 has no E8.75)
SEEDS = h28.SEEDS
OWN, CLEAN, TRUTHS, WORLDS = h28.OWN, h28.CLEAN, h28.TRUTHS, h28.WORLDS
FLOOR = harness.FLOOR
STRENGTH = 1.25  # OTS125's
COPY = "copy_D_f"  # the gate: row H29b's D_f (its organ:f), a reference
T10, T15 = h34.arm(1.0), h34.arm(1.5)  # the gate: row H34's two arms on D_f, references
GATE_BUILT = (COPY, T10, T15)
CONTROL, X, ANTI = "C", "X", "anti_X"
PLACEBOS = tuple(f"P{i}" for i in range(1, 10))
PLACEBO_SEEDS = dict(zip(PLACEBOS, range(20261251, 20261260), strict=True))
PERM_BINS = 10
COUNTED = (X, ANTI, *PLACEBOS)  # eleven tries
BUILT = (CONTROL, *COUNTED)
BAR, MIN_SEEDS = 0.30, 9
SHAPE_TOL = 1e-9
# a mean library size more than 1% from C's is VOID. The challenge wrote "median": at MERFISH depths (56 to 83 counts a
# cell) the integer median moves by one count, 1.2% to 2.4%, from count noise alone in many draws, which would void
# the read by construction; the mean (what keep_depth holds) moves about 0.1%. The median is kept as a diagnostic.
DEPTH_TOL = 0.01
SHAPE_TERMS, DE_TERMS = h28.SHAPE_TERMS, h28.DE_TERMS
H29B_READING_SHA256 = h34.H29B_READING_SHA256
H34_PARTS = HERE / "experiments" / "harness" / "results" / "h34_heart_organ_trend" / "parts" / "e95"
H34_PART_SHA256 = {  # row H34's committed parts on this transition (main at 397b0e1)
    "0-1": "f5b7f0023593ea95ffc6124497cf9428e0d3beee4444987502a9c41ea9591986",
    "2-3": "1f0f9dc9e8243241cf5e6f775a9d5924e29d2def487b9db1d28f8f6235692e9e",
    "4-5": "acf92531cfcd356320401e1a125eba3a74e48291385b4b790d2948edcb9165aa",
    "6-7": "c6cbc6dac3bc229a04b0e4e771c0e140904d4ec8776b52a81b9ecd09d39369bd",
    "8-9": "87b895886ef9b6b76b96d0cd813ff504b11bdb35e9c9b8dacdfdc2fbb3a86fb0",
}
G1_WORLDS = h34.G1_WORLDS  # own, CLEAN and the four single truths
G1_ARMS = h34.G1_ARMS
STAT_FIELDS = h34.STAT_FIELDS
ENTER, STOP, PASS, FAIL = "ENTER", "STOP", "PASS", "FAIL"
VOID, NOT_PASSED, PASSED = "VOID", "NOT PASSED", "PASSED AT PRACTICE"
STEP_EXIT = {ENTER: 0, PASS: 0, STOP: 3, FAIL: 3}
DIAGNOSTIC = "A labelled diagnostic on one truth: a pass is a priced question to the captain, never adoption."
NOTHING = f"What follows: nothing is built; no rerun under this row; {FILE} is not built."


# ------------------------------------------------------------------------------------------------ statements
def stop_line(where: str, why) -> str:
    return f"Try {ROW}: STOPPED AT {where} ({'; '.join(why)}): nothing further runs; {FILE} is not built."


def void_line(why) -> str:
    return f"Try {ROW}: VOID: {'; '.join(why)}; nothing is built."


def not_passed_line(bars) -> str:
    return f"Try {ROW}: NOT PASSED ({', '.join(bars)}): {FILE} is not built."


def t_stopped_line(why) -> str:
    return f"Try {ROW}: T STOPPED ({', '.join(why)}): {FILE} is not built."


FOLLOWS = {
    PASSED: (f"What follows: the outcome is committed into row {ROW}'s outcome cell word for word, then the test "
             f"profile step (T) runs at that commit; nothing is built before T passes."),
    NOT_PASSED: NOTHING,
    VOID: NOTHING,
}


# ------------------------------------------------------------------------------------------------ refusals
def git(*a) -> str:
    return subprocess.run(["git", "-C", str(HERE), *a], capture_output=True, text=True).stdout.strip()


def git_state() -> dict:
    return h29b.run_provenance()


def vec_home() -> Path:
    return Path(harness.__file__).resolve().parents[2]


def refuse_state(commit: str | None, what: str) -> None:
    """Every step's refusals, before anything loads: row H36 missing or malformed, vec from another checkout, tracked
    files that differ from the commit, a HEAD other than the commit given."""
    try:
        why = h36_row.row_reasons(TRIES.read_text(encoding="utf-8"))
    except OSError:
        why = ["TRIES.md cannot be read"]
    if why:
        sys.exit(f"{ROW} is refused ({what}): {'; '.join(why)} (a row comes before its read)")
    if vec_home() != HERE:
        sys.exit(f"{ROW} is refused ({what}): vec is imported from {vec_home()}, not from this checkout")
    state = git_state()
    if state.get("tree_clean") is not True:
        sys.exit(f"{ROW} is refused ({what}): tracked files differ from the commit")
    if not commit or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused ({what}): HEAD is not the commit {commit!r}")


def load_json(path: Path) -> dict | None:
    try:
        r = json.loads(Path(path).read_text())
        return r if isinstance(r, dict) else None
    except (OSError, ValueError):
        return None


def verdict_reasons(rec: dict | None, step: str, want: str, commit: str | None) -> list[str]:
    """The record of an earlier step: of this row, this step, verdict `want`, made at `commit` on a clean tree."""
    if not isinstance(rec, dict) or rec.get("row") != ROW or rec.get("step") != step:
        return [f"no {step} record"]
    if rec.get("verdict") != want:
        return [f"the {step} step said {rec.get('verdict')!r}, not {want}"]
    prov = rec.get("provenance") or {}
    if not commit or rec.get("commit") != commit or prov.get("commit") != commit or prov.get("tree_clean") is not True:
        return [f"the {step} record is not of the run commit"]
    return []


def _must_not_exist(*paths: Path) -> None:
    for p in paths:
        if Path(p).exists():
            sys.exit(f"{ROW} is refused: {p} exists (no step twice)")


def refuse(step: str, commit: str | None, out: Path | None = None, before: Path | None = None) -> None:
    """One step's refusals. `before`: the folder of the step it follows (entry for profile, profile for gate, gate for
    run); the read's own is the gate."""
    refuse_state(commit, step)
    want = {"profile": ("entry", ENTER), "gate": ("profile", PASS), "run": ("gate", PASS), "read": ("gate", PASS)}
    if step in want:
        name, verdict = want[step]
        why = verdict_reasons(load_json(Path(before) / f"{name}.json"), name, verdict, commit)
        if why:
            sys.exit(f"{ROW} is refused ({step}): {'; '.join(why)}")
    if out is not None:
        files = {"entry": ("entry.json",), "profile": ("profile.json", "practice_profile.json"),
                 "gate": ("gate.json",), "run": (PART,), "read": ("reading.json", "reading.txt")}.get(step, ())
        _must_not_exist(*(Path(out) / f for f in files))


# ------------------------------------------------------------------------------------------------ E and Q
def _record(step: str, commit: str, verdict: str, **kw) -> dict:
    return {"row": ROW, "axis": AXIS, "step": step, "commit": commit, "provenance": git_state(), "verdict": verdict,
            "source": SOURCE, **kw}


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "x", encoding="utf-8") as f:
        f.write(text)


def run_entry(out: Path, data_dir: Path, commit: str, read=None) -> int:
    refuse("entry", commit, out)
    genes = xpc.panel_genes()
    raw = (read or xpc.read_files)(data_dir, xpc.STAGES, genes)
    ent = xpc.entry(xpc.groups(raw, xpc.STAGES))
    verdict = ENTER if ent["ok"] else STOP
    _write(Path(out) / "entry.json", json.dumps(_record("entry", commit, verdict, entry=ent,
                                                       label="not a try: a look at E10.5's labels (H29b's caveat)"),
                                                indent=1))
    if verdict == STOP:
        print(stop_line("ENTRY", [f"under {xpc.MIN_CELLS} cells: {', '.join(ent['short'])}"]), flush=True)
    print(f"ENTRY {verdict}", flush=True)
    return STEP_EXIT[verdict]


def run_profile(out: Path, data_dir: Path, commit: str, entry_dir: Path, read=None) -> int:
    refuse("profile", commit, out, entry_dir)
    genes = xpc.panel_genes()
    stages = xpc.KIND_STAGES[xpc.PRACTICE]  # E8.5 and E9.5: no E10.5 file is read
    raw = (read or xpc.read_files)(data_dir, stages, genes)
    res = xpc.measure(xpc.groups(raw, stages), xpc.PRACTICE)
    files = {inst: {s: xpc.FILES[inst][s] for s in stages} for inst in xpc.INSTRUMENTS}
    text = xpc.dumps({**xpc.record(res, genes, files), "row": ROW, "commit": commit})
    _write(Path(out) / "practice_profile.json", text)
    verdict = PASS if res["passed"] else STOP
    _write(Path(out) / "profile.json", json.dumps(_record(
        "profile", commit, verdict, failed=res["failed"], practice_profile_sha256=xpc.sha256_text(text),
        label="one T1 reading (rows N59 to N61): E8.5 and E9.5 only, no E10.5 file"), indent=1))
    if verdict == STOP:
        print(stop_line("THE PRACTICE PROFILE", res["failed"]), flush=True)
    print(f"PROFILE {verdict}", flush=True)
    return STEP_EXIT[verdict]


def practice_delta(profile_dir: Path, genes) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """The committed practice profile on `genes` (by name): delta_p, its eligible genes and the reference-stage level
    the placebos are binned on, and the file's sha256 (the one the profile step recorded; refused otherwise)."""
    rec = load_json(Path(profile_dir) / "profile.json") or {}
    sha = rec.get("practice_profile_sha256")
    path = Path(profile_dir) / "practice_profile.json"
    prof = json.loads(path.read_bytes())
    d = xpc.load_profile(path, sha, prof.get("genes"), kind=xpc.PRACTICE)
    idx = {g: i for i, g in enumerate(prof["genes"])}
    if sorted(idx) != sorted(genes):
        raise ValueError("the practice profile does not name the transition's panel")
    order = np.array([idx[g] for g in genes])
    e = np.asarray(prof["eligible"], dtype=bool)[order]
    level = np.asarray(prof["level"], dtype=np.float64)[order]
    return d[order], e, level, sha


# ------------------------------------------------------------------------------------------------ the arms
def placebo_perm(level: np.ndarray, eligible: np.ndarray, seed: int) -> np.ndarray:
    """A permutation of the genes that moves every eligible gene to another eligible gene of the same level bin (ten
    bins of the eligible genes by rank of `level`, row H28's derangement) and keeps every other gene in place."""
    perm = np.arange(len(level))
    idx = np.flatnonzero(eligible)
    if idx.size:
        sub = h28.derangement(np.asarray(level, dtype=np.float64)[idx], f"{ROW}-{seed}", bins=PERM_BINS, seed=seed)
        perm[idx] = idx[sub]
    return perm


def arm_deltas(delta_p: np.ndarray, eligible: np.ndarray, level: np.ndarray) -> tuple[dict, dict]:
    """Each counted arm's profile: X delta_p, anti-X -delta_p, P_i delta_p deranged (placebo_perm); and each placebo's
    fixed points among the eligible genes."""
    out = {X: delta_p.copy(), ANTI: -delta_p}
    fixed = {}
    for name, seed in PLACEBO_SEEDS.items():
        p = placebo_perm(level, eligible, seed)
        out[name] = delta_p[p]
        fixed[name] = int((p[eligible] == np.flatnonzero(eligible)).sum())
    return out, fixed


def build_arms(ctx: harness.BuildContext, target: dict, deltas: dict) -> tuple[dict, dict]:
    """One seed: C (h28.trend_expression at 1.25 on D_f's context) and each arm, C's counts through
    vec.t2_build.apply_run_profile with the arm's profile at stream [build seed, 37]; each then D_f's size and half;
    and the identity checks the read's void conditions hold the arms to."""
    st_c, counts_c, rows_c, rec_c = h28.trend_expression(ctx, STRENGTH)
    pool, rows, _, coords = harness.draw_copy(ctx)
    stream = [ctx.seed, t2_build.RUN_PROFILE_STREAM]
    pre = {CONTROL: (st_c, counts_c)}
    for name, d in deltas.items():
        new, Xa, _ = t2_build.apply_run_profile(counts_c, d, stream)
        pre[name] = (harness.finish(ctx, pool, Xa, coords), new)
    # the arms' rendering at a zero profile is C itself, bit for bit (no count drawn, the same X and coordinates)
    zero = harness.finish(ctx, pool, t2_build.apply_run_profile(counts_c, np.zeros(counts_c.shape[1]), stream)[1],
                          coords)
    lib = {k: np.asarray(v[1].sum(1), dtype=np.float64).ravel() for k, v in pre.items()}
    checks = {"zero_profile_is_c": bool(np.array_equal(zero.X, st_c.X) and np.array_equal(zero.coords, st_c.coords)),
              "same_draw": bool(np.array_equal(rows, rows_c)),
              "pre_size_coords_equal": bool(all(np.array_equal(v[0].coords, st_c.coords) for v in pre.values())),
              # the depth guard gates on the mean (what keep_depth holds); the integer median is a diagnostic only
              "mean_library_ratio": {k: float(v.mean() / lib[CONTROL].mean()) for k, v in lib.items()},
              "median_library_ratio": {k: float(np.median(v) / np.median(lib[CONTROL])) for k, v in lib.items()}}
    out, recs = {}, {}
    for name, (st, _) in pre.items():
        out[name], recs[name] = h34.size_and_half(ctx, st, target, name)
    checks["post_size_coords_equal"] = bool(all(np.array_equal(v.coords, out[CONTROL].coords) for v in out.values()))
    checks["n_cells"] = {k: int(v.n) for k, v in out.items()}
    return out, {"checks": checks, "size": recs, "trend": rec_c, "draw": h34.draw_record(ctx, st_c)}


# ------------------------------------------------------------------------------------------------ the gate
def gate_seed(tr, s: int, data_dir: Path) -> dict:
    t0 = time.time()
    ctx = h34.context(tr, s)
    with harness.guarded_reads(harness.REPO_ROOT, data_dir):
        copy, rec, draw = h34.build_control(ctx, tr.target)
        preds = {COPY: copy}
        for strength, name in ((1.0, T10), (1.5, T15)):
            st, _, _, _ = h28.trend_expression(ctx, strength)
            preds[name], _ = h34.size_and_half(ctx, st, tr.target, name)
    scores, skills, groups = h34.score_arms(tr, s, preds)
    return {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "size_record": rec, "draw": draw,
            "scores": scores, "skills": skills, "groups": groups, "seconds": round(time.time() - t0, 1)}


def h34_committed_rows(root: Path = H34_PARTS) -> list[dict] | None:
    """Row H34's committed parts on this transition, merged by seed: only the files of the sha256 pinned here."""
    parts = []
    try:
        for seeds, sha in H34_PART_SHA256.items():
            raw = (Path(root) / seeds / h34.PART).read_bytes()
            if hashlib.sha256(raw).hexdigest() != sha:
                return None
            parts.append(json.loads(raw))
        return h26.merge_rows(parts)[0]
    except (OSError, ValueError, KeyError, TypeError):
        return None


def world_means(rows: list[dict], name: str) -> dict[str, float]:
    """The plain ten-seed means of `name` in the worlds own and CLEAN on this transition (h28.series)."""
    return {w: float(h28.series({TASK: rows}, (TASK,), ts, name).mean()) for w, ts in WORLDS}


def gate_reasons(rows: list[dict], facts: dict, committed: dict | None, h27_rows: list[dict] | None,
                 h34_rows: list[dict] | None) -> tuple[list[str], dict]:
    """Why the gate fails (empty when it passes) and what it found. G1 on this transition, exactly: (R), (D), (M) and
    (P) of row H34's G1 with this row's copy; G2: row H34's two arms' world means against its committed parts."""
    found: dict = {}
    why = []
    if [r["seed"] for r in rows] != list(SEEDS):
        return ["G1: the gate did not build every seed"], found
    if not h34.committed_ok(committed):
        why.append("G1: row H29b's reading.json is not at hand as the file its row records")
    else:
        found["R"] = h34.record_differs(TASK, facts, committed)
        why += [f"G1 (R): {', '.join(found['R'])} not row H29b's committed"] if found["R"] else []
        found["D"] = h34.pool_draw(rows) == committed["pools_short"][TASK][h34.POOL]
        why += [] if found["D"] else ["G1 (D): the copy's draw is not row H29b's committed"]
        unit = h29b.UNITS[TASK]
        mine = {w: float(h29b.series({TASK: rows}, unit, ts, COPY).mean()) for w, ts in G1_WORLDS}
        off = [w for w, _ in G1_WORLDS if mine[w] != committed["tables"][w][TASK]["mean"][h34.D_F]]
        found["M"] = {"means": len(G1_WORLDS), "unequal": off}
        why += [f"G1 (M): {len(off)} of {len(G1_WORLDS)} means are not row H29b's committed"] if off else []
        if h27_rows is None:
            why.append("G1 (P): row H27's committed parts are not at hand")
        else:
            off = []
            for other in G1_ARMS:
                for w, ts in G1_WORLDS:
                    got = h29b._stats(h29b.series({TASK: rows}, unit, ts, COPY)
                                      - h29b.series({TASK: h27_rows}, unit, ts, other))
                    want = committed["tables"][w][TASK]["stated"][h29b.pname(h34.D_F, other)]
                    off += [f"{w}/{other}/{k}" for k in STAT_FIELDS if got[k] != want[k]]
            found["P"] = {"numbers": len(G1_ARMS) * len(G1_WORLDS) * len(STAT_FIELDS), "unequal": off}
            why += [f"G1 (P): {len(off)} paired statistics are not row H29b's committed"] if off else []
    if h34_rows is None:
        why.append("G2: row H34's committed parts are not at hand as the files pinned here")
    else:
        off = [f"{name} {w}" for name in (T10, T15) for w, m in world_means(rows, name).items()
               if m != world_means(h34_rows, name)[w]]
        found["G2"] = {"means": 4, "unequal": off}
        why += [f"G2: {', '.join(off)} not row H34's committed mean"] if off else []
    return why, found


def run_gate(out: Path, data_dir: Path, commit: str, profile_dir: Path) -> int:
    refuse("gate", commit, out, profile_dir)
    committed, h27_all = h34.load_h29b_reading(), h34.load_h27_rows()
    h27_rows = None if h27_all is None else h27_all.get(TASK)
    tr = h34.load_transition(TASK, data_dir)
    rows = []
    for s in SEEDS:
        rows.append(gate_seed(tr, s, data_dir))
        print(f"{TASK} seed {s} ({rows[-1]['seconds']} s) done", flush=True)
    state = git_state()
    if state.get("tree_clean") is not True or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: tracked files or HEAD changed while the gate ran; no gate record is written")
    try:
        why, found = gate_reasons(rows, tr.facts, committed, h27_rows, h34_committed_rows())
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:
        why, found = [f"G1: a committed record cannot be read ({type(e).__name__}: {e})"], {}
    verdict = FAIL if why else PASS
    _write(Path(out) / "gate.json", json.dumps(_record(
        "gate", commit, verdict, why=why, found=found, h29b_reading_sha256=H29B_READING_SHA256,
        h34_part_sha256=H34_PART_SHA256, transition=TASK, facts=tr.facts, rows=rows,
        label="counting nothing: references only (row H29b's D_f, row H34's trend_s1.0 and trend_s1.5)"),
        indent=1, default=float))
    print(f"GATE {verdict}", flush=True)
    if why:
        print(void_line([f"the gate fails ({'; '.join(why)})"]), flush=True)
        print(f"{DIAGNOSTIC} {FOLLOWS[VOID]}", flush=True)
    return STEP_EXIT[verdict]


# ------------------------------------------------------------------------------------------------ the run
def run_part(seeds, out: Path, data_dir: Path, commit: str, gate_dir: Path, profile_dir: Path) -> dict:
    refuse("run", commit, out, gate_dir)
    tr = h34.load_transition(TASK, data_dir)
    genes = list(tr.cboard.genes)
    delta_p, eligible, level, sha = practice_delta(profile_dir, genes)  # read before the guard: results are refused
    deltas, fixed = arm_deltas(delta_p, eligible, level)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        ctx = h34.context(tr, s)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            preds, extra = build_arms(ctx, tr.target, deltas)
        scores, skills, groups = h34.score_arms(tr, s, preds)
        rows_out.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "scores": scores, "skills": skills,
                         "groups": groups, **extra, "seconds": round(time.time() - t0, 1)})
        print(f"{TASK} seed {s} ({rows_out[-1]['seconds']} s) done", flush=True)
    res = {"row": ROW, "axis": AXIS, "task": TASK, "arms": list(BUILT), "truths": list(TRUTHS), "seeds": list(seeds),
           "practice_profile_sha256": sha,
           "profile_checks": {"x_is_delta_p": bool(np.array_equal(deltas[X], delta_p)),
                              "anti_is_minus_delta_p": bool(np.array_equal(deltas[ANTI], -delta_p)),
                              "placebo_fixed_points": fixed},
           "pool_cells": int(tr.inputs["ref"].n_obs), "rows": rows_out, "source": SOURCE, "provenance": git_state()}
    path = Path(out) / PART
    _write(path, json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


def load_parts(root: Path) -> tuple[list[dict], list[str]]:
    parts, unreadable = [], []
    for p in sorted(Path(root).rglob(PART)):
        r = load_json(p)
        if r is None or r.get("axis") != AXIS or r.get("task") != TASK or not isinstance(r.get("rows"), list):
            unreadable.append(str(p))
        else:
            parts.append(r)
    return parts, unreadable


def _finite(v) -> bool:
    return isinstance(v, int | float) and not isinstance(v, bool) and math.isfinite(v)


def void_reasons(parts: list[dict], gate: dict | None, profile: dict | None, commit: str | None,
                 unreadable=()) -> list[str]:
    """Why the read is void, by the conditions code can hold: the gate's or the profile's record missing, failed or
    of another commit; a seed, truth or arm missing or doubled, a score not finite, a part unreadable; parts from
    more than one commit, another commit or an unclean tree; an arm that is not the arm defined (another draw,
    coordinates or cell number than C's, a shape term more than 1e-9 from C's, X's profile not delta_p or anti-X's
    not -delta_p, a placebo with a fixed point, a mean library more than 1% from C's, a profile not of the sha256
    the profile step recorded)."""
    why = [f"v: {w}" for w in verdict_reasons(gate, "gate", PASS, commit)]
    why += [f"v: {w}" for w in verdict_reasons(profile, "profile", PASS, commit)]
    why += [f"v: a part cannot be read: {p}" for p in unreadable]
    if not parts:
        return why + ["v: no part"]
    rows, twice = h26.merge_rows(parts)
    if twice:
        why.append(f"v: seeds {twice} occur in more than one part")
    if [r["seed"] for r in rows] != list(SEEDS):
        why.append(f"v: seeds {[r['seed'] for r in rows]}, not {list(SEEDS)}")
    if len({str(p["provenance"].get("commit")) for p in parts}) > 1:
        why.append("v: the parts come from more than one commit")
    if not all(p["provenance"].get("tree_clean") is True for p in parts):
        why.append("v: a part comes from an unclean tree")
    if not commit or any(str(p["provenance"].get("commit")) != commit for p in parts):
        why.append(f"v: a part comes from a commit other than the run commit {commit!r}")
    want_sha = (profile or {}).get("practice_profile_sha256")
    for p in parts:
        pc = p.get("profile_checks") or {}
        if p.get("practice_profile_sha256") != want_sha:
            why.append("v: a part's practice profile is not the file of the sha256 the profile step recorded")
        if pc.get("x_is_delta_p") is not True or pc.get("anti_is_minus_delta_p") is not True:
            why.append("v: X's profile is not delta_p or anti-X's is not -delta_p")
        if set(pc.get("placebo_fixed_points") or {}) != set(PLACEBOS) or any(pc["placebo_fixed_points"].values()):
            why.append("v: a placebo leaves a gene on its own value (or is missing)")
    for r in rows:
        tag = f"v: seed {r['seed']}"
        for field in ("scores", "skills"):
            gone = sorted({f"{q}/{a}" for q in TRUTHS for a in (FLOOR, *BUILT) if a not in r.get(field, {}).get(q, {})})
            if gone:
                why.append(f"{tag}: missing from {field}: {gone[:6]}")
        if any(not _finite(r["scores"][q][a]) for q in TRUTHS for a in (FLOOR, *BUILT)
               if a in r.get("scores", {}).get(q, {})):
            why.append(f"{tag}: a score is not a finite number")
        c = r.get("checks") or {}
        for k in ("zero_profile_is_c", "same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
            if c.get(k) is not True:
                why.append(f"{tag}: {k.replace('_', ' ')} fails")
        n = c.get("n_cells") or {}
        if set(n) != set(BUILT) or len(set(n.values())) != 1:
            why.append(f"{tag}: the arms' cell numbers are not C's")
        lib = c.get("mean_library_ratio") or {}
        if set(lib) != set(BUILT) or any(not (_finite(v) and abs(v - 1.0) <= DEPTH_TOL) for v in lib.values()):
            why.append(f"{tag}: a mean library size is more than 1% from C's")
        for q in TRUTHS:
            for name in COUNTED:
                for m in SHAPE_TERMS:
                    try:
                        ok = abs(r["skills"][q][name][m] - r["skills"][q][CONTROL][m]) <= SHAPE_TOL
                    except (KeyError, TypeError):
                        ok = False
                    if not ok:
                        why.append(f"{tag} {q}: {name}'s {m} is not C's")
    return why


# ------------------------------------------------------------------------------------------------ the bars
def series(rows: list[dict], truths, name: str, what: str = "score") -> np.ndarray:
    return h28.series({TASK: rows}, (TASK,), truths, name, what)


def placebo_mean(rows: list[dict], truths, what: str = "score") -> np.ndarray:
    return np.mean([series(rows, truths, p, what) for p in PLACEBOS], axis=0)


def bars(rows: list[dict]) -> dict:
    """B1 to B4 in each world (own, CLEAN), on this transition, every difference paired by seed."""
    out = {}
    for w, ts in WORLDS:
        x, c, a = series(rows, ts, X), series(rows, ts, CONTROL), series(rows, ts, ANTI)
        pm = placebo_mean(rows, ts)
        means = {p: float(series(rows, ts, p).mean()) for p in PLACEBOS}
        spread = float(np.std(list(means.values()), ddof=1))
        need = max(BAR, 2.0 * spread)
        b1, b2, b3 = x - c, x - pm, a - pm
        de = {m: float((series(rows, ts, X, m) - placebo_mean(rows, ts, m)).mean()) for m in DE_TERMS}
        out[w] = {
            "B1": {"mean": float(b1.mean()), "above": int((b1 > 0).sum()),
                   "holds": bool(b1.mean() >= BAR and (b1 > 0).sum() >= MIN_SEEDS)},
            "B2": {"mean": float(b2.mean()), "above": int((b2 > 0).sum()), "need": need, "placebo_sd": spread,
                   "x_mean": float(x.mean()), "placebo_means": means,
                   "holds": bool(b2.mean() >= need and (b2 > 0).sum() >= MIN_SEEDS
                                 and all(float(x.mean()) > m for m in means.values()))},
            "B3": {"mean": float(b3.mean()), "below": int((b3 < 0).sum()),
                   "holds": bool(b3.mean() <= -BAR and (b3 < 0).sum() >= MIN_SEEDS)},
            "B4": {**de, "holds": bool(all(v > 0 for v in de.values()))},
        }
    return out


def reading(parts: list[dict], gate: dict | None, profile: dict | None, commit: str | None, unreadable=()) -> dict:
    try:
        why = void_reasons(parts, gate, profile, commit, unreadable)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:
        why = [f"v: a part cannot be read ({type(e).__name__}: {e})"]
    res = {"row": ROW, "axis": AXIS, "commit": commit, "void": bool(why), "why": why, "diagnostic": DIAGNOSTIC}
    if why:
        res.update(outcome=VOID, line=void_line(why))
    else:
        rows = h26.merge_rows(parts)[0]
        b = bars(rows)
        failed = [f"{k} in {w}" for w, _ in WORLDS for k in ("B1", "B2", "B3", "B4") if not b[w][k]["holds"]]
        res.update(bars=b, failed=failed, outcome=NOT_PASSED if failed else PASSED,
                   line=not_passed_line(failed) if failed else h36_row.PASSED_AT_PRACTICE)
    res["follows"] = f"{DIAGNOSTIC} {FOLLOWS[res['outcome']]}"
    return res


def text(res: dict) -> str:
    lines = [res["line"], res["follows"],
             f"Row {ROW} (commit {res['commit']}): the XPC practice read on {TASK}, worlds own and CLEAN, seeds 0 "
             "to 9, "
             "released stages and the organisers' RNA at E8.5 and E9.5 only. No figure here is to be read beside a "
             "test score."]
    if res["void"]:
        return "\n".join([*lines, "VOID, the conditions found:", *res["why"]]) + "\n"
    for w, _ in WORLDS:
        for k, v in res["bars"][w].items():
            lines.append(f"{w} {k}: {'holds' if v['holds'] else 'does not hold'} "
                         + json.dumps({a: b for a, b in v.items() if a != 'holds'}, default=float))
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ T
def refuse_test(commit: str, reading_dir: Path, out: Path | None = None, profile_path: Path = PROFILE) -> dict:
    """T's refusals: the state; the read's reading tracked and unchanged at HEAD, PASSED AT PRACTICE; the row's
    outcome cell holding its statement word for word; no test profile yet. Returns the reading."""
    out = Path(reading_dir) / "test_profile" if out is None else Path(out)
    refuse_state(commit, "test-profile")
    path = Path(reading_dir) / "reading.json"
    rel = path.resolve().relative_to(HERE) if path.resolve().is_relative_to(HERE) else None
    if rel is None or not git("ls-files", "--", str(rel)):
        sys.exit(f"{ROW} is refused (test-profile): the read's reading is not committed ({path})")
    res = load_json(path)
    if not res or res.get("row") != ROW or res.get("void") is not False or res.get("outcome") != PASSED:
        sys.exit(f"{ROW} is refused (test-profile): the read did not say {PASSED}")
    cell = h36_row.outcome_cell(TRIES.read_text(encoding="utf-8")) or ""
    if h36_row.PASSED_AT_PRACTICE not in cell or res.get("line") != h36_row.PASSED_AT_PRACTICE:
        sys.exit(f"{ROW} is refused (test-profile): row {ROW}'s outcome cell does not hold the read's statement")
    if any(e in cell for e in h36_row.ENDINGS):
        sys.exit(f"{ROW} is refused (test-profile): row {ROW}'s outcome cell holds a statement that ends the line")
    _must_not_exist(Path(profile_path), out)
    return res


def run_test(data_dir: Path, commit: str, reading_dir: Path, out: Path | None = None, profile_path: Path = PROFILE,
             read=None) -> int:
    out = Path(reading_dir) / "test_profile" if out is None else Path(out)
    refuse_test(commit, reading_dir, out, profile_path)
    genes = xpc.panel_genes()
    raw = (read or xpc.read_files)(data_dir, xpc.STAGES, genes)
    gs = xpc.groups(raw, xpc.STAGES, (xpc.DISSECTED,))
    alt = xpc.groups(raw, xpc.STAGES, xpc.ALL_E105_LIBRARIES)
    res = xpc.measure(gs, xpc.TEST, alt)
    files = {inst: dict(xpc.FILES[inst]) for inst in xpc.INSTRUMENTS}
    body = xpc.dumps({**xpc.record(res, genes, files), "row": ROW, "commit": commit})
    sha = xpc.sha256_text(body)
    verdict = PASS if res["passed"] else STOP
    _write(Path(profile_path) if verdict == PASS else Path(out) / "xpc_profile_not_passed.json", body)
    line = h36_row.pass_statement(sha) if verdict == PASS else t_stopped_line(res["failed"])
    _write(Path(out) / "test_profile.json", json.dumps(_record(
        "test-profile", commit, verdict, failed=res["failed"], profile_sha256=sha, line=line,
        label="one heart reading (E10.5) and one T1 reading (library 3), after the read's outcome was committed"),
        indent=1))
    print(line, flush=True)
    return STEP_EXIT[verdict]


# ------------------------------------------------------------------------------------------------ command line
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    steps = {}
    for name in ("entry", "profile", "gate", "run", "read", "test-profile", "refusals"):
        p = sub.add_parser(name)
        p.add_argument("--commit", required=True, help="the commit HEAD must be (full hash)")
        p.add_argument("--data", default=str(harness.DATA_DIR))
        p.add_argument("--out", default=str(RESULTS))
        steps[name] = p
    steps["run"].add_argument("--seeds", required=True)
    steps["refusals"].add_argument("--step", required=True,
                                   choices=("entry", "profile", "gate", "run", "read", "test-profile"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if a.cmd == "refusals":
        step = a.step
        if step == "test-profile":
            refuse_test(a.commit, out)
        elif step == "run":
            refuse(step, a.commit, None, out / "gate")
        else:
            sub_out = {"entry": out / "entry", "profile": out / "profile", "gate": out / "gate", "read": out}[step]
            before = {"profile": out / "entry", "gate": out / "profile", "read": out / "gate"}.get(step)
            refuse(step, a.commit, sub_out, before)
        return
    if a.cmd == "entry":
        sys.exit(run_entry(out / "entry", Path(a.data), a.commit))
    if a.cmd == "profile":
        sys.exit(run_profile(out / "profile", Path(a.data), a.commit, out / "entry"))
    if a.cmd == "gate":
        sys.exit(run_gate(out / "gate", Path(a.data), a.commit, out / "profile"))
    if a.cmd == "run":
        run_part(h28.payoff.parse_seeds(a.seeds), out / "parts" / a.seeds, Path(a.data), a.commit, out / "gate",
                 out / "profile")
        return
    if a.cmd == "read":
        refuse("read", a.commit, out, out / "gate")
        parts, unreadable = load_parts(out / "parts")
        res = reading(parts, load_json(out / "gate" / "gate.json"), load_json(out / "profile" / "profile.json"),
                      a.commit, unreadable)
        _write(out / "reading.json", json.dumps(res, indent=1, default=float))
        _write(out / "reading.txt", text(res))
        print(res["line"])
        print(res["follows"])
        return
    sys.exit(run_test(Path(a.data), a.commit, out))


if __name__ == "__main__":
    main()
