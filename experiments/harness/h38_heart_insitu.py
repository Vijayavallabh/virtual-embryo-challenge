#!/usr/bin/env python
"""Row H38 of experiments/harness/TRIES.md: the counted read of heart file OIC125 (row R70, the captain's order of
2026-11-25; HEART_MOON1_design.md P1 with every fix text of HEART_MOON1_challenge.md section 5). Ten counted heart
tries on one truth (E10.5): ISC, anti-ISC and three placebos on the two backtests T2h/<=E9.5->E10.5 and
T2h/<=E8.75->E10.5, merged seed by seed into fold E10.5; worlds own and CLEAN (row H26's); seeds 0 to 9. A labelled
diagnostic under rule 9 J5. It builds no test file and changes no file, rule or default.

The arms, each on row H35's D_f context (D's pool the reference input), then D_f's size and half (h34.size_and_half):
  C         OTK125's practice form: the category form at 1.25 (vec.trend_forms.form_step), row H35's category arm,
            rebuilt by the gate to equal H35's committed parts exactly (a re-read: not counted)
  ISC       C with each category's practice chord P added where P exists (vec.trend_forms.insitu_step; FIX 1)
  anti-ISC  C with -P
  PL1-PL3   C with P deranged among its genes within ten level bins of the category's reference level (log mean CP10k
            + 0.5 over D's pool's cells of the category), seeds 20261126 to 20261128, no fixed point
The practice chords are G0's (experiments/heart_ext/oic_g0.py): chords_E9.5.json for the first backtest and
chords_E8.75.json for the second, each the file of the sha256 that row H38's outcome cell records in G0's PASS
statement.

  gate      counting nothing: per transition and seed, C built and scored; G1 holds C's and the floor's scores equal to
            row H35's committed category arm and floor in every seed and truth, and the in-situ step with its chord out
            equal to the category form bit for bit; the J5 floor (twice the larger world's paired seed SE of H35's
            category minus its control on fold E10.5); exit 0 PASS, 3 FAIL (VOID before any counted arm)
  run       one transition, a seed range: the six arms built and scored, FIX 6's faults, libraries, draws; one part
  read      the void conditions first, then the bars, on fold E10.5 in both worlds:
              S1 ISC - C >= max(+0.30, J5 floor), at least 9 of 10 seeds above 0
              S2 ISC - each placebo >= max(+0.20, J5 floor), at least 8 of 10 seeds above 0
              S3 anti-ISC - C <= -0.10
              S4 ISC - C >= 0 on each transition
              S5 (replacing the design's S5, which could bite only through the shape group's nbhd gate, the shape
                 terms being C's) the mean skill differences ISC - C in de_score and de_direction both above 0
            PASS when every bar holds in both worlds; STOP (IDLE) when |ISC - C| < 0.10 on the fold in both worlds;
            STOP otherwise; VOID decided first. Writes reading.json and reading.txt; prints the outcome line and what
            follows.
  j5-floor  prints the J5 floor from row H35's committed parts (practice figures only)
  refusals  a step's refusals and nothing else
Every step is refused before anything loads unless TRIES.md holds exactly one row H38 of eleven cells whose first cell
is h38_row.FIRST_CELL and whose outcome cell holds G0's PASS statement and no read outcome, the chords are committed
files of the sha256 it records, the tracked tree is clean, HEAD is the run commit and vec is this checkout's own.
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

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "src"))

from vec import harness, t2_build, trend_forms  # noqa: E402
from vec.core import counts_of  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h35 = _load("h35_heart_trend_forms", HERE / "experiments" / "harness" / "h35_heart_trend_forms.py")
h34, h28, h26, payoff = h35.h34, h35.h28, h35.h26, h35.payoff
hrow = _load("h38_row", HERE / "experiments" / "harness" / "h38_row.py")
oic = _load("oic_g0", HERE / "experiments" / "heart_ext" / "oic_g0.py")
oc = oic.oc

ROW, RULING_ROW, FILE = hrow.ROW, hrow.RULING, hrow.FILE
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h38-heart-insitu"
GATE_AXIS = AXIS + ":gate"
PART = "h38.json"
GATE_JSON, GATE_TXT = "gate.json", "gate.txt"
RESULTS = HERE / "experiments" / "harness" / "results" / "h38_heart_insitu"
CHORDS = HERE / "experiments" / "heart_ext" / "results"
SOURCE = "experiments/harness/h38_heart_insitu.py (row H38)"

STRENGTH = 1.25
TASK_KEYS = h35.TASK_KEYS
REF_OF = {"T2h/<=E9.5->E10.5": "E9.5", "T2h/<=E8.75->E10.5": "E8.75"}
SEEDS, TRUTHS, WORLDS, UNITS, FOLD = h35.SEEDS, h35.TRUTHS, h35.WORLDS, h35.UNITS, h35.FOLD
TOL = h35.TOL
SHAPE_TERMS, DE_TERMS, TERMS = h35.SHAPE_TERMS, h35.DE_TERMS, h35.TERMS
FLOOR = harness.FLOOR
C, ISC, ANTI = "C", "ISC", "anti-ISC"
PLACEBOS = ("PL1", "PL2", "PL3")
PL_SEEDS = (20261126, 20261127, 20261128)
COUNTED = (ISC, ANTI, *PLACEBOS)
BUILT = (C, *COUNTED)
H35_ARM = "category"  # row H35's name of C's arm
S1_BAR, S1_SEEDS = 0.30, 9
S2_BAR, S2_SEEDS = 0.20, 8
S3_BAR = -0.10
IDLE_BAR = 0.10
LIB_TOL = 0.01
J5_MULT = 2.0
PASS, FAIL, STOP, IDLE, VOID = "PASS", "FAIL", "STOP", "STOP (IDLE)", "VOID"
GATE_EXIT = {PASS: 0, FAIL: 3}


def pname(a: str, b: str) -> str:
    return f"{a} - {b}"


# ---- refusals -----------------------------------------------------------------------------------------------------
def git_state() -> dict:
    return h35.git_state()


def vec_home() -> Path:
    return Path(harness.__file__).resolve().parents[2]


def tracked_clean(path: Path) -> bool:
    rel = str(Path(path).resolve().relative_to(HERE))
    t = subprocess.run(["git", "-C", str(HERE), "ls-files", "--error-unmatch", rel], capture_output=True)
    d = subprocess.run(["git", "-C", str(HERE), "status", "--porcelain", "--", rel], capture_output=True, text=True)
    return t.returncode == 0 and not d.stdout.strip()


def chord_shas(text: str | None = None) -> tuple[dict | None, list[str]]:
    """The chords' sha256 by reference from row H38's G0 statement, each a committed file of that sha256."""
    if text is None:
        try:
            text = TRIES.read_text(encoding="utf-8")
        except OSError:
            text = ""
    shas, why = hrow.g0_shas(text)
    if shas is None:
        return None, why
    out = dict(zip(("E10.5", "E9.5", "E8.75"), shas, strict=True))
    for ref, sha in out.items():
        p = CHORDS / f"chords_{ref}.json"
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != sha:
            why.append(f"{p.name} is not the file of row {ROW}'s sha256")
        elif not tracked_clean(p):
            why.append(f"{p.name} is not committed as it stands")
    return (out if not why else None), why


def refuse_state(commit: str | None, what: str = "run") -> dict:
    shas, why = chord_shas()
    if why:
        sys.exit(f"{ROW} is refused before its {what}: " + "; ".join(why))
    if vec_home() != HERE:
        sys.exit(f"{ROW} is refused: vec is imported from {vec_home()}, not from this checkout")
    state = git_state()
    if state.get("tree_clean") is not True:
        sys.exit(f"{ROW} is refused: tracked files differ from the commit")
    if not commit or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit {commit!r}")
    return shas


def _load_json(path: Path) -> dict | None:
    try:
        r = json.loads(Path(path).read_text())
        return r if isinstance(r, dict) else None
    except (OSError, ValueError):
        return None


def _of_commit(rec: dict, commit: str | None) -> bool:
    prov = rec.get("provenance") or {}
    return bool(commit and rec.get("commit") == commit and prov.get("commit") == commit
                and prov.get("tree_clean") is True)


def gate_record_reasons(gate: dict | None, commit: str | None) -> list[str]:
    if not isinstance(gate, dict) or gate.get("axis") != GATE_AXIS:
        return ["v3: no gate record"]
    if gate.get("verdict") != PASS:
        return [f"v3: the gate said {gate.get('verdict')!r}, not {PASS}"]
    if not _of_commit(gate, commit):
        return ["v3: the gate record is not of the run commit"]
    return []


def refuse_gate(out: Path, commit: str | None) -> None:
    refuse_state(commit, "gate")
    if (Path(out) / GATE_JSON).exists() or (Path(out) / GATE_TXT).exists():
        sys.exit(f"{ROW} is refused: {out} already holds a gate record")


def refuse_run(path: Path | None, commit: str | None, gate_dir: Path) -> None:
    refuse_state(commit)
    if gate_record_reasons(_load_json(Path(gate_dir) / GATE_JSON), commit):
        sys.exit(f"{ROW} is refused: the gate of {commit!r} has not said {PASS}")
    if path is not None and Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


def refuse_read(out: Path, commit: str | None) -> None:
    refuse_state(commit, "read")
    if (Path(out) / "reading.json").exists() or (Path(out) / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")


# ---- the arms -----------------------------------------------------------------------------------------------------
def ref_levels(ctx) -> dict[str, np.ndarray]:
    """Per category, log(mean CP10k + 0.5) over D's pool's cells of that category (the placebos' level basis)."""
    pool = ctx.inputs[ctx.task.reference_role]
    genes = list(pool.var_names)
    cats = t2_build.t2_chamber_categories(np.asarray(pool.obs[t2_build.T2_LABEL_KEY]).astype(str))
    C = counts_of(pool, genes)
    out = {}
    for k in oc.CATS:
        m = cats == k
        out[k] = np.log(oc.cp10k(C[np.flatnonzero(m)]).mean(0) + oc.PSEUDO) if m.any() else np.zeros(len(genes))
    return out


def variants(P: dict[str, np.ndarray], levels: dict[str, np.ndarray]) -> tuple[dict, dict]:
    """ISC's chord, its negative and the three placebos; and each placebo's fixed points (must be none)."""
    out = {ISC: {k: v.copy() for k, v in P.items()}, ANTI: {k: -v for k, v in P.items()}}
    fixed = {}
    for name, seed in zip(PLACEBOS, PL_SEEDS, strict=True):
        d, fx = {}, 0
        for k, v in P.items():
            idx = np.flatnonzero(np.isfinite(v))
            perm = oc.derange_within_bins(levels[k][idx], key=k, seed=seed)
            fx += int((perm == np.arange(len(idx))).sum())
            w = np.full(len(v), np.nan)
            w[idx] = v[idx][perm]
            d[k] = w
        out[name], fixed[name] = d, fx
    return out, fixed


def lib_mean(counts) -> float:
    return float(np.asarray(counts.sum(1)).ravel().mean())


def build_arms(ctx, target: dict, P: dict[str, np.ndarray]) -> tuple[dict, dict]:
    """Every arm of one seed, through D_f's size and half, and the checks the read's v5 holds them to."""
    pre = {C: trend_forms.form_step(ctx, STRENGTH, "category")}
    var, fixed = variants(P, ref_levels(ctx))
    base = trend_forms.form_arguments(ctx, STRENGTH, "category")
    faults = {}
    for name in COUNTED:
        args = trend_forms.insitu_arguments(ctx, STRENGTH, var[name])
        faults[name] = trend_forms.insitu_checks(ctx, args, base, var[name])
        pre[name] = trend_forms.insitu_step(ctx, STRENGTH, var[name], check=False)
    rows = pre[C][2]
    checks = {"same_draw": bool(all(np.array_equal(v[2], rows) for v in pre.values())),
              "pre_size_coords_equal": bool(all(np.array_equal(v[0].coords, pre[C][0].coords) for v in pre.values())),
              "lib_mean": {k: lib_mean(v[1]) for k, v in pre.items()}, "faults": faults, "fixed_points": fixed}
    out, recs = {}, {}
    for name, (st, _, _, rec) in pre.items():
        out[name], srec = h34.size_and_half(ctx, st, target, name)
        recs[name] = {"size": srec, "carried": rec.get("carried", []), "form": rec.get("form"),
                      "types": rec.get("types"), "strength": rec.get("strength")}
    checks["post_size_coords_equal"] = bool(all(np.array_equal(v.coords, out[C].coords) for v in out.values()))
    checks["n_cells"] = {k: int(v.n) for k, v in out.items()}
    return out, {"records": recs, "checks": checks}


def transition_chords(key: str, shas: dict, genes) -> dict[str, np.ndarray]:
    ref = REF_OF[key]
    return oic.load_chords(CHORDS / f"chords_{ref}.json", sha256=shas[ref], genes=list(genes))


def run_transition(key: str, seeds, out: Path, data_dir: Path, commit: str | None, gate_dir: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse_run(path, commit, gate_dir)
    shas = refuse_state(commit)
    if key not in TASK_KEYS:
        sys.exit(f"{key}: not a transition of row {ROW}")
    tr = h35.load_transition(key, data_dir)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        ctx = h35.context(tr, s)
        P = transition_chords(key, shas, ctx.inputs[ctx.task.reference_role].var_names)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            preds, extra = build_arms(ctx, tr.target, P)
        scores, skills, groups = h34.score_arms(tr.base, s, preds)
        rows_out.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, **extra, "scores": scores,
                         "skills": skills, "groups": groups, "seconds": round(time.time() - t0, 1)})
        print(f"{key} seed {s} ({rows_out[-1]['seconds']} s) done", flush=True)
    res = {"label": f"row {ROW}: C (OTK125's practice form) and the in-situ arms on row H29b's D_f form; a labelled "
                    "diagnostic on one truth", "axis": AXIS, "task": tr.key, "ref": REF_OF[key],
           "chords_sha256": shas[REF_OF[key]], "arms": list(BUILT), "seeds": list(seeds), "rows": rows_out,
           "strength": STRENGTH, "source": SOURCE, "provenance": git_state()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


# ---- the gate -----------------------------------------------------------------------------------------------------
def h35_rows() -> dict[str, list[dict]] | None:
    parts, bad = h35.load_parts(h35.RESULTS)
    if bad or set(parts) != set(TASK_KEYS):
        return None
    return {k: h26.merge_rows(parts[k])[0] for k in TASK_KEYS}


def j5_floor(rows35: dict[str, list[dict]] | None) -> float | None:
    """Twice the larger world's paired seed SE of row H35's category minus its control on fold E10.5."""
    if rows35 is None:
        return None
    se = [h34.compare(rows35, TASK_KEYS, ts, H35_ARM, h35.CONTROL)["score"]["se"] for _, ts in WORLDS]
    return float(J5_MULT * max(se))


def gate_seed(tr, s: int, data_dir: Path, h35_by_seed: dict | None) -> dict:
    ctx = h35.context(tr, s)
    with harness.guarded_reads(harness.REPO_ROOT, data_dir):
        st, _, _, _ = trend_forms.form_step(ctx, STRENGTH, "category")
        z = trend_forms.insitu_step(ctx, STRENGTH, None, zero=True, check=False)[0]
        pred, _ = h34.size_and_half(ctx, st, tr.target, C)
    scores, _, _ = h34.score_arms(tr.base, s, {C: pred})
    ident = bool(np.array_equal(z.X, st.X) and np.array_equal(z.coords, st.coords))
    g1 = None
    if h35_by_seed is not None and s in h35_by_seed:
        theirs = h35_by_seed[s]["scores"]
        g1 = bool(all(scores[q][C] == theirs[q][H35_ARM] and scores[q][FLOOR] == theirs[q][FLOOR] for q in TRUTHS))
    return {"seed": s, "scores": scores, "zero_identity": ident, "g1": g1}


def run_gate(out: Path, data_dir: Path, commit: str | None) -> int:
    out = Path(out)
    refuse_gate(out, commit)
    rows35 = h35_rows()
    rows, why = {}, []
    if rows35 is None:
        why.append("G1: row H35's committed parts are missing or unreadable")
    for key in TASK_KEYS:
        tr = h35.load_transition(key, data_dir)
        by = {r["seed"]: r for r in rows35[key]} if rows35 else None
        rows[key] = [gate_seed(tr, s, data_dir, by) for s in SEEDS]
        for r in rows[key]:
            if r["g1"] is not True:
                why.append(f"G1: {key} seed {r['seed']}: C or the floor is not row H35's committed arm")
            if not r["zero_identity"]:
                why.append(f"{key} seed {r['seed']}: the in-situ step with its chord out is not the category form")
    floor = j5_floor(rows35)
    verdict = FAIL if why else PASS
    rec = {"row": ROW, "axis": GATE_AXIS, "verdict": verdict, "why": why, "j5_floor": floor,
           "bars": {"S1": max(S1_BAR, floor or 0.0), "S2": max(S2_BAR, floor or 0.0)}, "rows": rows,
           "commit": commit, "provenance": git_state()}
    out.mkdir(parents=True, exist_ok=True)
    h35._write_new(out / GATE_JSON, json.dumps(rec, indent=1, default=float))
    h35._write_new(out / GATE_TXT, f"GATE {ROW}: {verdict}." + ("" if not why else " " + "; ".join(why)) + "\n")
    print(f"GATE {ROW}: {verdict}.", flush=True)
    return GATE_EXIT[verdict]


# ---- the read -----------------------------------------------------------------------------------------------------
def load_parts(root: Path) -> tuple[dict[str, list[dict]], list[str]]:
    out: dict[str, list[dict]] = {}
    bad = []
    for p in sorted(Path(root).rglob(PART)):
        try:
            r = json.loads(p.read_text())
            if not (isinstance(r, dict) and r.get("axis") == AXIS and isinstance(r.get("rows"), list)):
                raise ValueError
        except (OSError, ValueError):
            bad.append(str(p))
            continue
        out.setdefault(r["task"], []).append(r)
    return out, bad


def _finite(v) -> bool:
    return isinstance(v, int | float) and not isinstance(v, bool) and math.isfinite(v)


def void_reasons(parts: dict, gate: dict | None, commit: str | None, shas: dict | None, unreadable=()) -> list[str]:
    why = gate_record_reasons(gate, commit) + [f"v4: a part cannot be read: {p}" for p in unreadable]
    if shas is None:
        why.append("v3: row H38's G0 statement or its chords are missing")
    every = [p for ps in parts.values() for p in ps]
    complete = True
    for key in TASK_KEYS:
        ps = parts.get(key)
        if not ps:
            why.append(f"v4: missing transition {key}")
            complete = False
            continue
        rows, twice = h26.merge_rows(ps)
        if twice:
            why.append(f"v4: {key}: seeds {twice} occur in more than one part")
            complete = False
        if [r["seed"] for r in rows] != list(SEEDS):
            why.append(f"v4: {key}: seeds {[r['seed'] for r in rows]}, not {list(SEEDS)}")
            complete = False
        for r in rows:
            gone = sorted({f"{q}/{a}" for q in TRUTHS for a in (FLOOR, *BUILT)
                           if not _finite(r.get("scores", {}).get(q, {}).get(a))})
            if gone:
                why.append(f"v4: {key} seed {r['seed']}: missing or not finite: {gone[:6]}")
                complete = False
        if shas is not None and any(p.get("chords_sha256") != shas[REF_OF[key]] for p in ps):
            why.append(f"v3: {key}: a part was made with another chord file")
    if set(parts) - set(TASK_KEYS):
        why.append(f"v4: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if every:
        if len({p["provenance"]["commit"] for p in every}) > 1:
            why.append("v4: the parts come from more than one commit")
        if not all(p["provenance"]["tree_clean"] is True for p in every):
            why.append("v4: a part comes from an unclean tree")
        if not commit or any(str(p["provenance"]["commit"]) != commit for p in every):
            why.append(f"v4: a part comes from a commit other than the run commit {commit!r}")
    if not complete or why:
        return why
    gate_scores = {k: {r["seed"]: r["scores"] for r in gate["rows"][k]} for k in TASK_KEYS}
    for key in TASK_KEYS:
        for r in h26.merge_rows(parts[key])[0]:
            tag = f"v5: {key} seed {r['seed']}"
            c = r["checks"]
            for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
                if c.get(k) is not True:
                    why.append(f"{tag}: {k.replace('_', ' ')} fails")
            if set(c["n_cells"]) != set(BUILT) or len(set(c["n_cells"].values())) != 1:
                why.append(f"{tag}: the arms' cell numbers are not C's")
            for name in COUNTED:
                if not abs(c["lib_mean"][name] / c["lib_mean"][C] - 1.0) <= LIB_TOL:
                    why.append(f"{tag}: {name}'s mean library size is more than 1% from C's")
                for f in c["faults"].get(name, []):
                    why.append(f"{tag}: {name}: {f}")
                for q in TRUTHS:
                    for m in SHAPE_TERMS:
                        if not abs(r["skills"][q][name][m] - r["skills"][q][C][m]) <= TOL:
                            why.append(f"{tag} {q}: {name}'s {m} is not C's")
            for name, n in c["fixed_points"].items():
                if n:
                    why.append(f"{tag}: {name} has {n} fixed points")
            for q in TRUTHS:
                for a in (C, FLOOR):
                    if r["scores"][q][a] != gate_scores[key][r["seed"]][q][a]:
                        why.append(f"{tag} {q}: the run step's {a} is not the gate's")
    return why


def compare(rows, keys, truths, a, b) -> dict:
    return h34.compare(rows, keys, truths, a, b)


PAIRS = ((ISC, C), (ANTI, C), *((ISC, p) for p in PLACEBOS), *((p, C) for p in PLACEBOS), (C, FLOOR), (ISC, FLOOR))


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: {pname(a, b): compare(rows, ks, ts, a, b) for a, b in PAIRS} for u, ks in UNITS.items()}
            for w, ts in WORLDS}


def bars(tab: dict, s1: float, s2: float) -> dict:
    """The bars in each world on fold E10.5 (S4 on each transition)."""
    out = {}
    for w, _ in WORLDS:
        f = tab[w][FOLD]
        d = f[pname(ISC, C)]["score"]
        out[w] = {"S1": bool(d["mean"] >= s1 and d["ahead"] >= S1_SEEDS),
                  "S2": bool(all(f[pname(ISC, p)]["score"]["mean"] >= s2
                                 and f[pname(ISC, p)]["score"]["ahead"] >= S2_SEEDS for p in PLACEBOS)),
                  "S3": bool(f[pname(ANTI, C)]["score"]["mean"] <= S3_BAR),
                  "S4": bool(all(tab[w][k][pname(ISC, C)]["score"]["mean"] >= 0 for k in TASK_KEYS)),
                  "S5": bool(all(f[pname(ISC, C)][m] > 0 for m in DE_TERMS)),
                  "idle": bool(abs(d["mean"]) < IDLE_BAR)}
    return out


FOLLOWS = {PASS: ("What follows: one priced question to the captain with both worlds' readings as ranges and the row's "
                  "ninth cell's caveats; on his word experiments/builds/heart_oic_build.sh builds OIC125 twice under "
                  "ruling row R70 after OTK125's identity build; never adoption."),
           STOP: ("What follows: nothing is built under row R70; the line ends there, with no second read under this "
                  "or another row."),
           IDLE: ("What follows: nothing is built under row R70 (the chords change nothing that the read can see); "
                  "the line ends there."),
           VOID: ("What follows: nothing is built under row R70 and no rerun under this row; the main session states "
                  "the void conditions.")}


def reading(parts: dict, gate: dict | None, commit: str | None, shas: dict | None, unreadable=()) -> dict:
    try:
        why = void_reasons(parts, gate, commit, shas, unreadable)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError, ZeroDivisionError) as e:
        why = [f"v4: a part cannot be read ({type(e).__name__}: {e})"]
    res = {"row": ROW, "commit": commit, "void": bool(why), "why": why}
    if why:
        outcome = VOID
    else:
        rows = {k: h26.merge_rows(parts[k])[0] for k in TASK_KEYS}
        tab = tables(rows)
        s1, s2 = gate["bars"]["S1"], gate["bars"]["S2"]
        b = bars(tab, s1, s2)
        passed = all(b[w][s] for w, _ in WORLDS for s in ("S1", "S2", "S3", "S4", "S5"))
        idle = all(b[w]["idle"] for w, _ in WORLDS)
        outcome = PASS if passed else IDLE if idle else STOP
        res.update(bars=b, tables=tab, s1_bar=s1, s2_bar=s2, j5_floor=gate.get("j5_floor"))
    res.update(outcome=outcome, line=hrow.outcome_statement(outcome), follows=FOLLOWS[outcome])
    return res


def text(res: dict) -> str:
    lines = [res["line"], res["follows"],
             f"Row {ROW} (commit {res['commit']}): the counted read of OIC125's practice form on row H29b's D_f, one "
             "truth (E10.5), worlds own and CLEAN, seeds 0 to 9; practice figures, never beside a test score."]
    if res["void"]:
        return "\n".join([*lines, "VOID, the conditions found:", *res["why"]]) + "\n"
    lines.append(f"Bars: S1 {res['s1_bar']:.3f}, S2 {res['s2_bar']:.3f} (J5 floor {res['j5_floor']}).")
    for w, _ in WORLDS:
        lines.append(f"{w}: " + ", ".join(f"{k} {v}" for k, v in res["bars"][w].items()))
        for a, b in PAIRS:
            s = res["tables"][w][FOLD][pname(a, b)]["score"]
            lines.append(f"  {pname(a, b)}: {s['mean']:+.3f} [{s['lo']:+.3f}, {s['hi']:+.3f}] {s['ahead']}/10 above "
                         f"(se {s['se']:.3f})")
    return "\n".join(lines) + "\n"


# ---- commands -----------------------------------------------------------------------------------------------------
def cmd_gate(a) -> None:
    sys.exit(run_gate(Path(a.out), Path(a.data), a.commit))


def cmd_run(a) -> None:
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS}")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data), a.commit, Path(a.gate))


def cmd_read(a) -> None:
    out = Path(a.out)
    refuse_read(out, a.commit)
    shas, _ = chord_shas()
    parts, bad = load_parts(Path(a.runs))
    res = reading(parts, _load_json(Path(a.gate) / GATE_JSON), a.commit, shas, bad)
    out.mkdir(parents=True, exist_ok=True)
    h35._write_new(out / "reading.json", json.dumps(res, indent=1, default=float))
    h35._write_new(out / "reading.txt", text(res))
    print(res["line"])
    print(res["follows"])


def cmd_floor(a) -> None:
    f = j5_floor(h35_rows())
    if f is None:
        sys.exit("row H35's committed parts are missing or unreadable")
    print(f"J5 floor (twice the larger world's paired seed SE of H35's category - control on fold E10.5): {f:.4f}")


def cmd_refusals(a) -> None:
    if a.step == "run":
        if a.gate is None:
            sys.exit("refusals --step run needs --gate")
        refuse_run(None, a.commit, Path(a.gate))
    elif a.out is None:
        sys.exit(f"refusals --step {a.step} needs --out")
    elif a.step == "gate":
        refuse_gate(Path(a.out), a.commit)
    else:
        refuse_read(Path(a.out), a.commit)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--commit", required=True)
    g.add_argument("--data", default=str(harness.DATA_DIR))
    r = sub.add_parser("run")
    r.add_argument("--task", required=True)
    r.add_argument("--seeds", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--gate", required=True)
    r.add_argument("--commit", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--gate", required=True)
    d.add_argument("--out", default=str(RESULTS))
    d.add_argument("--commit", required=True)
    sub.add_parser("j5-floor")
    f = sub.add_parser("refusals")
    f.add_argument("--step", required=True, choices=("gate", "run", "read"))
    f.add_argument("--commit", required=True)
    f.add_argument("--out")
    f.add_argument("--gate")
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "read": cmd_read, "j5-floor": cmd_floor, "refusals": cmd_refusals}[a.cmd](a)


if __name__ == "__main__":
    main()
