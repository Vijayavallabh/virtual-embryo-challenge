#!/usr/bin/env python
"""Try K31 (experiments/harness/TRIES.md, committed before any arm is scored): G1 of row R72's sign cube, the counted
read on T3's practice tasks (T3_MOON5_design.md 2.6 as fixes F2, F6, F7, F9 and F10 word it). A labelled diagnostic:
it asks only whether the stage and the Mab21l2 axes are large enough, at the test's reference and stage, that their
unknown signs move the test-like truth a board point beyond a random direction of their size. The sign of a
difference never picks a vertex; the Gata4 axis is unreadable here (the truth is the Gata4 run).

The analogue cube without a Gata4 table: A' = (S_M + Y_M) / 2, X_S' = (S_M - Y_M) / 2, X_M' = l_M - S_M, so the base
A' + X_S' + X_M' = l_M is K18's `barred_form` factor exactly. Plain world only (K8's harness and scorer, the reference
as released), K18's stream `t3_s2_form` and cells, seeds 0-9, tasks g1 and g2 (WT E8.75 to Gata4 replicates 1 and 2)
merged as one truth (J4: each seed's difference is the mean of the two tasks'). Every arm goes through K18's
barred_form path (build_submission.prediction_inputs with the barred s1's spec, construct, assemble) with its own
factor in place of phi_M, so every arm holds the same drawn rows. Per seed, B's columns are m_L * phi_M for every
label with a cell in that seed's draw (F1: "G1 uses phi_lM and K18's cells"; every label, as row R72 holds the cube's
files), Q projects off B, and every arm
is held as a cube vertex is (sign_cube.flip_factor: k for the SD of ln phi_M, the clip, the den constant over WT
E8.75's mean within-cell proportions, the Newton hold of each label's expected library):
  barred_form   phi_M itself: the guard arm, not counted while inside its guard
  cube_np       the analogue vertex np, the stage axis reflected: base - (1 + k) Q(X_S') + (k - 1) Q(X_M'), k solved
                for SD(ln phi_M) (at k = 1 it is base - 2 Q(X_S')), held (try)
  cube_pn       the analogue vertex pn, the Mab21l2 axis reflected: base + (k - 1) Q(X_S') - (1 + k) Q(X_M'), k
                solved for SD(ln phi_M) (at k = 1 it is base - 2 Q(X_M')), held (try)
  null_S        cube_np with Q(perm_s(Q(X_S'))) in place of Q(X_S'), k solved and held as cube_np is (at k = 1 it is
                base - 2 Q(perm_s(Q(X_S')))) (try; ten permutations, one counted arm)
  null_M        cube_pn with Q(perm_s(Q(X_M'))) in place of Q(X_M'), k solved and held as cube_pn is (try)
perm_s is vec.harness.t3_decile_permutation within deciles of the reference's mean expression (K18's basis) seeded
20261113 + s at draw seed s. An axis that G0 held at '+' is not read: its arm and null are not run and not counted.
An arm whose factor cannot be made at a seed (sign_cube.CubeRefused) is recorded as not made, and its axis cannot pass.

The reading (fixed in row K31 before any score):
  VOID        base off K18's committed plain barred_form (its task files at their pins) by > 0.1 in the ten-seed
              mean or > 0.5 in any seed, on either task (the base then counts one try); files of mixed commits or
              machines, a dirty tree, a seed file without every arm, G0's reading off its record
  PASS_j      |mean Delta_j| >= 1.0 board point (Delta_j,s = arm_j - base at seed s, the mean over g1 and g2),
              >= 8 of 10 seeds on the mean's side, and |mean Delta_j| above the second largest of the ten
              |Delta_null_j,s| (it exceeds at least nine of the ten; a tie fails)
  CONTINUE    iff some PASS_j: the cube is built over the passing axes plus X_G if G0 left it, each vertex only if G0
              leaves it; a failed axis is held at '+'
  STOP        otherwise (nothing changes: no axis moves the test-like truth a point beyond a random direction of its
              size); also when G0 said STOP; nothing is built
  INCOMPLETE  a seed file or run_started.json missing: not a reading (nothing opened beyond the file list)

  row                      prints the manifest sha256 row K31 records (2 when it holds none or several)
  dryrun --out FILE        J7: the guard applied to K21's committed reproduction of barred_form (reading.json's plain
                           guard distances against K18's committed rows); exit 0 holds, 6 fails
  start --out R             R/run_started.json, written once before the parts start (the launcher's run step), after
                           run's own checks: exit 0 written or already this run's, 2 refused
  run --task g1|g2 --seeds A-B --out R --data D --run-dir RUN --manifest M
                           one plain task, seed by seed, into R/plain_<task>/seed_<s>.json (needs R/g0.json, G0's
                           CONTINUE at this commit and machine; R/run_started.json is checked, and written once
                           through a per-process temporary name when it is absent)
  read --runs R --out OUT  the fixed reading into OUT/reading.json and OUT/reading.txt (and G0's reading copied there):
                           exit 0 a reading, 3 INCOMPLETE, 4 a reading already there (never overwritten), 2 refused,
                           1 a crash (masked)
No log holds a score: the harness's per-seed line is replaced by a progress line, and a crash's traceback has every
decimal number masked.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in (str(HERE), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import cube_tables as ct  # noqa: E402
import k7_wt_revert as k7  # noqa: E402  (the fold's reference and pool)
import k18_built_files as k18  # noqa: E402  (stream, tasks, barred spec and factor, stamps, folder rules)
import sign_cube as sc  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

ROW = "K31"
ROW_SCRIPT = "`experiments/t3_decomp/k31_cube_read.py`"
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"
MANIFEST_RE = re.compile(r"cube tables manifest sha256 `([0-9a-f]{64})`")
STREAM = k18.STREAM  # "t3_s2_form"
BASE = k18.BARRED  # "barred_form"
READ_AXES = ("S", "M")
ARM = {"S": "cube_np", "M": "cube_pn"}
NULL = {"S": "null_S", "M": "null_M"}
SIGNS = {"S": (-1, 1), "M": (1, -1)}  # in READ_AXES order
PERM_SEED = 20261113  # null_j at draw seed s: PERM_SEED + s
PASS_MEAN, PASS_SIDE = 1.0, 8
GUARD_MEAN, GUARD_SEED = k18.GUARD_MEAN, k18.GUARD_SEED  # 0.1 and 0.5
TASKS = k18.TASKS  # g1 and g2
K18_TASKS = HERE / "results" / "k18_built_files" / "tasks"
K18_SHA256 = {"plain_g1.json": "5ec09ef583862ab47ab078e72c7aed807c89b4d9c234c57e3b0727ff9e460ed9",
              "plain_g2.json": "cfc5ede7e7185ba5e16d4774d975b9e140e60bbd29736fb85c0fdb8e99f4634f"}
K21_READING = HERE / "results" / "k21_factor_forms" / "reading.json"
RESULTS = HERE / "results" / "k31_cube_read"
MANIFEST = ct.DEFAULT_OUT / ct.MANIFEST
LINE_VOID = "Try K31: VOID: {why}; nothing is built."
LINE_STOP = ("Try K31: STOP: no axis moves the test-like truth a point beyond a random direction of its size; nothing "
             "is built.")
LINE_STOP_G0 = "Try K31: STOP (G0): {why}; nothing is read or built."
LINE_GO = "Try K31: CONTINUE (PASS {passed}): the cube is built over {axes}: {vertices}."
INCOMPLETE = ("not a reading: a seed file or run_started.json is missing (a crash, a refusal such as a dirty tree, or "
              "a kill); nothing was read. Rerun the missing part at the run's commit on the same machine; after a code "
              "fix, rerun into a fresh folder")
EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_THERE, EXIT_GUARD = 2, 3, 4, 6

_BARRED: dict = {}  # the barred factor and spec, WT E8.75's genes (k18.load_barred_for_plain)
_CUBE: dict = {}  # the two tables and den in WT E8.75's gene order, and the axes read
_FOLD: dict = {}


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row and G0


def row_k31(path: Path | None = None) -> str:
    """Row K31's recorded manifest sha256: exactly one TRIES.md line starting '| K31 ' that names this script and
    holds exactly one 'cube tables manifest sha256 `<sha256>`'; Refused otherwise."""
    p = TRIES_MD if path is None else Path(path)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {ROW} ")] if p.is_file() \
        else []
    if len(lines) != 1 or ROW_SCRIPT not in lines[0]:
        raise ct.Refused(f"TRIES.md holds no single row {ROW} naming {ROW_SCRIPT}; the row is committed before any run")
    shas = set(MANIFEST_RE.findall(lines[0]))
    if len(shas) != 1:
        raise ct.Refused(f"row {ROW} does not record exactly one cube tables manifest sha256")
    return shas.pop()


def read_g0(path: Path, st: dict | None = None) -> dict:
    """G0's reading (cube_tables.py g0); with `st`, refused unless it was made at st's commit on st's machine from a
    clean tree."""
    p = Path(path)
    if not p.is_file():
        raise ct.Refused(f"{p}: no G0 reading")
    g0 = json.loads(p.read_text())
    prov = g0.get("provenance") or {}
    if st is not None and (prov.get("commit") != st["commit"] or prov.get("machine") != st["machine"]
                           or prov.get("tree_clean") is not True):
        raise ct.Refused("G0's reading is not of this commit and machine, or not from a clean tree")
    if g0.get("outcome") not in ("CONTINUE", "STOP"):
        raise ct.Refused("G0's reading has no outcome")
    return g0


def axes_read(g0: dict) -> tuple[str, ...]:
    return tuple(j for j in READ_AXES if j in (g0.get("axes_left") or []))


def arms_for(axes) -> tuple[str, ...]:
    return (BASE, *[a for j in axes for a in (ARM[j], NULL[j])])


def tries_for(axes) -> list[str]:
    return [a for j in axes for a in (ARM[j], NULL[j])]


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def load_cube(data: Path, run_dir: Path, decision: Path, manifest: Path, sha: str, axes) -> dict:
    """The barred factor and spec (k18.load_barred_for_plain), the tables S_M and Y_M at the manifest row K31 records,
    WT E8.75's den (build_outside_form.wt_den), all over WT E8.75's genes; held for the plain world's builds."""
    rec = k18.load_barred_for_plain(Path(data), Path(run_dir), Path(decision))
    _BARRED.clear()
    _BARRED.update(k18._BARRED)
    genes = list(_BARRED["genes"])
    tabs, man, man_sha = ct.load_tables(Path(manifest), genes, sha)
    _, den = ct._wt_den(Path(data) / ct.WT_FILES["E8.75"], genes)
    ax = sc.analogue_axes(np.log(_BARRED["phi"]), np.log(tabs["S_M"]), np.log(tabs["Y_M"]))
    _CUBE.clear()
    _CUBE.update(genes=genes, X={j: ax[j] for j in READ_AXES}, den=den, axes=tuple(axes), manifest_sha256=man_sha)
    _FOLD.clear()
    return rec


def fold(ctx) -> dict:
    """The fold's base, axes and den in the reference's gene order, and the permutation basis (K18's: the reference
    pool's mean expression), made once per prepared fold."""
    prep = k7._prepared(ctx)
    if _FOLD.get("prep") is not prep:
        if not _BARRED or not _CUBE:
            raise RuntimeError("the barred factor and the cube tables are not loaded (run loads them first)")
        ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
        genes = list(ref.genes)
        F = {"phi": k18.in_gene_order(_BARRED["phi"], _CUBE["genes"], genes),
             "X": {j: k18.in_gene_order(_CUBE["X"][j], _CUBE["genes"], genes) for j in READ_AXES},
             "den": k18.in_gene_order(_CUBE["den"], _CUBE["genes"], genes),
             "basis": harness.t3_mean_expression(pool), "spec": _BARRED["spec"]}
        _FOLD.clear()
        _FOLD.update(prep=prep, F=F, seeds={})
    return _FOLD["F"]


def null_axes(j: str, q: dict, Q, basis, seed: int) -> list[np.ndarray]:
    """null_j's projected axes: Q(perm_s(Q(X_j'))) in place of Q(X_j'), the other axis as the arm has it."""
    z = Q(harness.t3_decile_permutation(q[j], basis, PERM_SEED + int(seed)))
    return [z, q["M"]] if j == "S" else [q["S"], z]


def seed_factors(ref, n: int, seed: int, F: dict, axes) -> dict:
    """One seed's arms: the rows barred_form draws (construct with no factor on a fresh stable_rng(seed, STREAM)), the
    labels of them (every label with a drawn cell) and their rows m_L, B and Q, then every arm through
    sign_cube.flip_factor."""
    hit = _FOLD["seeds"].get(int(seed))
    if hit is not None:
        return hit
    bs = harness._t3_head()
    kw, _ = bs.prediction_inputs(ref, F["spec"], {}, None)
    drawn = bs.construct(ref, n, bs.t3d.stable_rng(seed, STREAM), **kw)
    if not np.all(drawn.source == "wt"):
        raise RuntimeError("barred_form's cells come from the reference pool only")
    held, M = sc.label_rows(ref.pool.C[drawn.row], drawn.lab)
    B = sc.basis(M, F["phi"])
    Q = sc.projector(B)
    q = {j: Q(F["X"][j]) for j in READ_AXES}
    phis, recs, not_made = {}, {}, {}
    for j in axes:
        for arm, qs in ((ARM[j], [q["S"], q["M"]]), (NULL[j], null_axes(j, q, Q, F["basis"], seed))):
            try:
                phis[arm], r = sc.flip_factor(F["phi"], qs, SIGNS[j], M, F["den"], B)
                recs[arm] = {k: r[k] for k in ("k", "sd_ln_factor", "max_abs_library_dev", "den_dev", "min", "max",
                                               "factor_sha256")} | {"rounds": r["hold"]["rounds"]}
            except sc.CubeRefused as e:
                not_made[arm] = str(e)
    out = {"phi": phis, "records": recs, "not_made": not_made, "held_labels": int(len(held)),
           "sd_q": {j: sc.sd(q[j]) for j in READ_AXES}}
    _FOLD["seeds"][int(seed)] = out
    return out


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = k7._prepared(ctx)
    F = fold(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    n = harness.count_for(ctx.board, pool.n_obs)
    phi = F["phi"]
    if arm != BASE:
        phi = seed_factors(ref, n, ctx.seed, F, _CUBE["axes"])["phi"].get(arm, F["phi"])  # not made: recorded
    kw, _ = bs.prediction_inputs(ref, F["spec"], {}, phi)
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    built = bs.construct(ref, n, bs.t3d.stable_rng(ctx.seed, STREAM), **kw)
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates(axes) -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params={"arm": arm, "stream": STREAM, "perm_seed": PERM_SEED if arm.startswith("null")
                                         else None, "row": ROW},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms_for(axes))]


# ---------------------------------------------------------------------------------------------------------------------
# Running


def progress_log(task: str):
    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"K31 plain {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


def parse_seeds(text_: str | None) -> list[int]:
    if text_ is None:
        return list(SEEDS)
    out: list[int] = []
    for part in str(text_).split(","):
        a, _, b = part.partition("-")
        lo, hi = int(a), int(b) if b else int(a)
        out += list(range(lo, hi + 1))
    if not out or any(s not in SEEDS for s in out) or len(set(out)) != len(out):
        raise ct.Refused(f"--seeds {text_!r}: distinct members of {list(SEEDS)}")
    return sorted(out)


def checked_axes(out: Path, st: dict) -> tuple[str, tuple[str, ...]]:
    """The run's checks before anything is scored: row K31's manifest sha256, G0's CONTINUE at this commit and machine
    on that manifest, and an axis left to read. Returns (the sha256, the axes); ct.Refused otherwise."""
    sha = row_k31()
    g0 = read_g0(Path(out) / "g0.json", st)
    if g0["outcome"] != "CONTINUE" or g0.get("manifest_sha256") != sha:
        raise ct.Refused("G0 said STOP, or read another manifest: nothing is scored")
    axes = axes_read(g0)
    if not axes:
        raise ct.Refused("G0 left neither S nor M: nothing is scored")
    return sha, axes


def started_record(st: dict, axes) -> dict:
    return {**st, "row": ROW, "axes": list(axes), "tries": tries_for(axes)}


def publish_once(path: Path, raw: bytes) -> bool:
    """`raw` written to `path` only when `path` does not exist, atomically and safely against parallel writers: a
    temporary file named for this process (its pid), then os.link to the final name, which fails when the name exists
    (exclusive creation where the file system has no hard links); this process's temporary name is then removed. True
    when this call wrote the file. Run k31a's four parts wrote one shared temporary name, and one part's rename
    removed it from under another."""
    path = Path(path)
    tmp = path.parent / f".{path.name}.{os.getpid()}.part"
    tmp.write_bytes(raw)
    try:
        os.link(tmp, path)
        return True
    except FileExistsError:
        return False
    except OSError:  # no hard links here: exclusive creation, then the bytes
        try:
            with open(path, "xb") as fh:
                fh.write(raw)
            return True
        except FileExistsError:
            return False
    finally:
        tmp.unlink(missing_ok=True)


def write_started(out: Path, st: dict, axes) -> dict:
    """run_started.json once per run folder (publish_once), then read back: ct.Refused unless it is of this commit and
    machine, from a clean tree, row K31 and these axes. Returns the record read."""
    path = Path(out) / "run_started.json"
    if not path.exists():
        publish_once(path, (json.dumps(started_record(st, axes), indent=1) + "\n").encode())
    rec = json.loads(path.read_text())
    if (rec.get("commit") != st["commit"] or rec.get("machine") != st["machine"] or rec.get("tree_clean") is not True
            or rec.get("row") != ROW or rec.get("axes") != list(axes)):
        raise ct.Refused("run_started.json is of another commit, machine or axes: a fresh run folder is needed")
    return rec


def cmd_start(a) -> int:
    """The run step's start record, written once by the launcher before the parts start (run_k31.sh): run's own checks
    (a clean tree, row K31, G0's CONTINUE at this commit and machine, an axis left), then write_started. Exit 0
    written or already this run's, 2 refused."""
    st = k18.stamp()
    if not st["tree_clean"]:
        print(k18.DIRTY, file=sys.stderr)
        return EXIT_REFUSED
    try:
        _, axes = checked_axes(Path(a.out), st)
        Path(a.out).mkdir(parents=True, exist_ok=True)
        rec = write_started(Path(a.out), st, axes)
    except ct.Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return EXIT_REFUSED
    print(f"run_started.json: commit {rec['commit']}, machine {rec['machine']}, axes {', '.join(rec['axes'])}",
          flush=True)
    return 0


def cmd_run(a) -> int:
    st = k18.stamp()
    if not st["tree_clean"]:
        print(k18.DIRTY, file=sys.stderr)
        return EXIT_REFUSED
    out = Path(a.out)
    try:
        sha, axes = checked_axes(out, st)
        seeds = parse_seeds(a.seeds)
        if a.task not in TASKS:
            raise ct.Refused("--task takes g1 or g2")
        folder = out / f"plain_{a.task}"
        folder.mkdir(parents=True, exist_ok=True)
        k18.check_folder(folder, st)
        write_started(out, st, axes)  # the launcher's start step wrote it; written here only when run alone
    except ct.Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return EXIT_REFUSED
    try:
        rec = load_cube(Path(a.data), Path(a.run_dir), Path(a.decision), Path(a.manifest), sha, axes)
        for s in seeds:
            path = folder / f"seed_{s}.json"
            if path.exists():
                print(f"K31 plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            res = harness.run_task(TASKS[a.task], candidates(axes), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log(a.task))
            got, now = res["provenance"], k18.stamp()
            if got.get("tree_clean") is not True or got.get("commit") != st["commit"] or not now["tree_clean"] \
                    or now["commit"] != st["commit"]:
                print("the tree changed while the part ran: no task file written (INCOMPLETE)", file=sys.stderr)
                return EXIT_REFUSED
            res["provenance"]["machine"] = st["machine"]
            sf = _FOLD["seeds"].get(int(s), {})
            res |= {"world": "plain", "task_source": "experiments/t3_decomp/k31_cube_read.py", "seed": s,
                    "axes": list(axes), "manifest_sha256": sha, "barred_factor": rec,
                    "factors": {k: v for k, v in sf.items() if k != "phi"}}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        return 1
    return 0


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def guard_distance(rows: list[dict], ref_rows: list[dict]) -> dict | None:
    """The base against K18's committed barred_form, seed by seed (k18.distance)."""
    return k18.distance(rows, BASE, ref_rows, BASE)


def pass_axis(delta, delta_null) -> dict:
    """PASS_j on ten per-seed differences (F2 and F6's bars)."""
    d, dn = np.asarray(delta, dtype=np.float64), np.abs(np.asarray(delta_null, dtype=np.float64))
    mean = float(d.mean())
    side = int(((np.sign(d) == np.sign(mean)) & (d != 0)).sum()) if mean != 0 else 0
    second = float(np.sort(dn)[-2]) if dn.size >= 2 else float("inf")
    above = bool(abs(mean) > second)
    return {"mean": mean, "seeds_on_side": side, "null_abs": [float(x) for x in dn],
            "null_second_largest": second, "above_null": above,
            "passes": bool(abs(mean) >= PASS_MEAN and side >= PASS_SIDE and above)}


def outcome_vertices(passed: list[str], g0: dict) -> tuple[list[str], list[str]]:
    """The axes the cube is built over (the passing ones, plus G if G0 left it) and the vertices built (those
    reflecting only those axes and left by G0)."""
    axes = [j for j in sc.AXES if j in passed or (j == "G" and "G" in (g0.get("axes_left") or []))]
    left = set(g0.get("vertices_left") or [])
    return axes, [v for v in sc.kept_vertices(axes) if v in left]


def reading(runs: Path, refs: Path = K18_TASKS) -> dict:
    runs = Path(runs)
    g0p = runs / "g0.json"
    if not g0p.is_file():
        return {"incomplete": True, "void": False, "missing": ["g0.json"], "outcome": INCOMPLETE}
    g0_raw = g0p.read_bytes()
    g0 = json.loads(g0_raw)
    base = {"g0": {"sha256": hashlib.sha256(g0_raw).hexdigest(), "outcome": g0.get("outcome"),
                   "line": g0.get("line"), "axes_left": g0.get("axes_left"), "vertices_left": g0.get("vertices_left"),
                   "provenance": g0.get("provenance")}}
    if g0.get("outcome") == "STOP":
        line = LINE_STOP_G0.format(why=str(g0.get("line", "")).rstrip("."))
        return base | {"incomplete": False, "void": False, "why": [], "outcome": "STOP", "line": line,
                       "tries_counted": [], "build": []}
    axes = axes_read(g0)
    files = ["run_started.json"] + [f"plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
    missing = [f for f in files if not (runs / f).exists()]
    if missing:
        return base | {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    why: list[str] = []
    started = json.loads((runs / "run_started.json").read_text())
    if started.get("axes") != list(axes):
        why.append("run_started.json's axes are not G0's")
    rows: dict = {n: [] for n in TASKS}
    prov, shas = {}, {}
    for n in TASKS:
        for s in SEEDS:
            f = f"plain_{n}/seed_{s}.json"
            raw = (runs / f).read_bytes()
            r = json.loads(raw)
            shas[f] = hashlib.sha256(raw).hexdigest()
            prov[f] = {k: (r.get("provenance") or {}).get(k) for k in ("commit", "tree_clean", "machine")}
            if r.get("seed") != s or r.get("task", {}).get("key") != TASKS[n].key:
                why.append(f"{f}: another seed or task")
                continue
            if r.get("axes") != list(axes) or r.get("manifest_sha256") != g0.get("manifest_sha256"):
                why.append(f"{f}: other axes or another manifest than G0's")
            seed_rows = r.get("rows") or []
            if len(seed_rows) != 1 or seed_rows[0].get("seed") != s:
                why.append(f"{f}: not one row of its seed")
                continue
            q = seed_rows[0]
            miss = [x for x in (harness.FLOOR, *arms_for(axes)) if x not in q.get("scores", {})]
            if miss:
                why.append(f"{f}: missing arms {miss}")
                continue
            nm = (r.get("factors") or {}).get("not_made") or {}
            rows[n].append(q | {"not_made": nm})
    for f, p in prov.items():
        if p["tree_clean"] is not True:
            why.append(f"{f}: a dirty tree")
    commits = {p["commit"] for p in prov.values()} | {started.get("commit")}
    machines = {p["machine"] for p in prov.values()} | {started.get("machine")}
    g0prov = g0.get("provenance") or {}
    if len(commits) > 1 or g0prov.get("commit") not in commits:
        why.append(f"results from different commits {sorted(map(str, commits | {g0prov.get('commit')}))}")
    if len(machines) > 1 or g0prov.get("machine") not in machines:
        why.append(f"results from different machines {sorted(map(str, machines | {g0prov.get('machine')}))}")
    guard = {}
    for n in TASKS:
        f = Path(refs) / f"plain_{n}.json"
        if not f.is_file() or ct.sha256(f) != K18_SHA256[f"plain_{n}.json"]:
            why.append(f"K18's committed plain_{n}.json is off its pin")
            continue
        if len(rows[n]) == len(SEEDS):
            guard[n] = guard_distance(sorted(rows[n], key=lambda x: x["seed"]), json.loads(f.read_text())["rows"])
            if not k18.guard_passes(guard[n]):
                d = guard[n]
                why.append(f"guard: plain {n}: {BASE} is off K18's committed rows" + (
                    f" by {d['mean_diff']:+.3f} in the mean and {d['max_abs_diff']:.3f} in one seed" if d else
                    " (a seed is missing)") + f" (limits {GUARD_MEAN} and {GUARD_SEED})")
    tries = tries_for(axes)
    out = base | {"incomplete": False, "files_sha256": shas, "provenance": prov, "guard": guard,
                  "axes_read": list(axes), "run_started": started}
    if why:
        guard_failed = any(w.startswith("guard:") for w in why)
        return out | {"void": True, "why": why, "outcome": "VOID", "line": LINE_VOID.format(why="; ".join(why)),
                      "tries_counted": tries + ([BASE] if guard_failed else []), "build": []}
    by_seed = {s: {n: next(x for x in rows[n] if x["seed"] == s) for n in TASKS} for s in SEEDS}
    axes_out, passed = {}, []
    for j in axes:
        made = all(ARM[j] not in by_seed[s][n]["not_made"] and NULL[j] not in by_seed[s][n]["not_made"]
                   for s in SEEDS for n in TASKS)
        d = [float(np.mean([by_seed[s][n]["scores"][ARM[j]] - by_seed[s][n]["scores"][BASE] for n in TASKS]))
             for s in SEEDS]
        dn = [float(np.mean([by_seed[s][n]["scores"][NULL[j]] - by_seed[s][n]["scores"][BASE] for n in TASKS]))
              for s in SEEDS]
        res = pass_axis(d, dn)
        res |= {"made": made, "delta": d, "delta_null": dn,
                "per_task": {n: {"arm_mean": float(np.mean([by_seed[s][n]["scores"][ARM[j]] -
                                                            by_seed[s][n]["scores"][BASE] for s in SEEDS])),
                                 "null_mean": float(np.mean([by_seed[s][n]["scores"][NULL[j]] -
                                                             by_seed[s][n]["scores"][BASE] for s in SEEDS]))}
                             for n in TASKS}}
        res["passes"] = bool(res["passes"] and made)
        axes_out[j] = res
        if res["passes"]:
            passed.append(j)
    if not passed:
        return out | {"void": False, "why": [], "axes": axes_out, "outcome": "STOP", "line": LINE_STOP,
                      "tries_counted": tries, "build": []}
    built_axes, verts = outcome_vertices(passed, g0)
    line = LINE_GO.format(passed=" and ".join(passed), axes=", ".join(built_axes), vertices=", ".join(verts) or
                          "no vertex (G0 left none of them)")
    return out | {"void": False, "why": [], "axes": axes_out, "outcome": "CONTINUE", "line": line, "passed": passed,
                  "built_axes": built_axes, "tries_counted": tries, "build": verts}


def text(res: dict) -> str:
    lines = [res.get("line") or res["outcome"]]
    if res.get("incomplete"):
        return "\n".join([*lines, "missing: " + ", ".join(res["missing"])]) + "\n"
    lines.append(f"G0: {res['g0'].get('line')} (g0.json sha256 {res['g0']['sha256']})")
    lines.append("tries counted: " + (", ".join(res.get("tries_counted") or []) or "none"))
    for n, d in (res.get("guard") or {}).items():
        lines.append(f"guard plain {n}: " + ("n/a" if d is None else
                                             f"mean {d['mean_diff']:+.4f}, largest seed {d['max_abs_diff']:.4f}"))
    for j, r in (res.get("axes") or {}).items():
        lines.append(f"axis {j}: mean Delta {r['mean']:+.3f}, {r['seeds_on_side']} of {len(SEEDS)} seeds on its side, "
                     f"second largest |null| {r['null_second_largest']:.3f}, made {r['made']}: "
                     f"{'PASS' if r['passes'] else 'no pass'}")
        lines.append("  per seed: " + ", ".join(f"{x:+.2f}" for x in r["delta"]))
        lines.append("  null per seed: " + ", ".join(f"{x:+.2f}" for x in r["delta_null"]))
    if res.get("void"):
        lines.append("VOID: " + "; ".join(res["why"]))
    lines.append("build: " + (", ".join(res.get("build") or []) or "nothing"))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> int:
    out = Path(a.out)
    if (out / "reading.json").exists() or (out / "reading.txt").exists():
        print(f"refused: {out} holds a reading already (never overwritten)", file=sys.stderr)
        return EXIT_THERE
    try:
        res = reading(Path(a.runs))
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        return 1
    if res.get("incomplete"):
        print(t, end="")
        return EXIT_INCOMPLETE
    out.mkdir(parents=True, exist_ok=True)
    for f in ("g0.json", "g0.txt"):
        if (Path(a.runs) / f).is_file():
            shutil.copyfile(Path(a.runs) / f, out / f)
    (out / "reading.json").write_text(js)
    (out / "reading.txt").write_text(t)
    print(t, end="")
    return 0


def dry_run(reading_path: Path = K21_READING, refs: Path = K18_TASKS) -> dict:
    """J7: this row's guard (k18.guard_passes, K17's tolerance) applied to K21's committed reproduction of barred_form
    against K18's committed plain rows (K21's reading.json, its plain guard distances), with K18's files at their
    pins. Reads committed results only; not a try."""
    r = json.loads(Path(reading_path).read_text())
    pins = {f: ct.sha256(Path(refs) / f) == h for f, h in K18_SHA256.items()}
    out = {}
    for n in TASKS:
        d = (r.get("guard") or {}).get(f"plain {n}")
        out[n] = {"distance": d, "passes": k18.guard_passes(d)}
    return {"k21_reading": str(reading_path), "k21_reading_sha256": ct.sha256(reading_path), "k18_pins": pins,
            "tasks": out, "passes": bool(all(pins.values()) and all(v["passes"] for v in out.values())),
            "could_fail_there": "no: K21 read distance 0 on every seed (the base arm is K18's own code path)"}


def cmd_dryrun(a) -> int:
    try:
        res = dry_run()
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        return 1
    Path(a.out).write_text(json.dumps(res, indent=1) + "\n")
    for n, v in res["tasks"].items():
        d = v["distance"]
        print(f"J7 dry run, plain {n}: " + ("n/a" if d is None else
                                            f"mean {d['mean_diff']:+.4f}, largest seed {d['max_abs_diff']:.4f}")
              + f": {'holds' if v['passes'] else 'fails'}")
    return 0 if res["passes"] else EXIT_GUARD


def cmd_row(a) -> int:
    """Prints the manifest sha256 row K31 records (exit 2 when the row does not hold exactly one)."""
    try:
        print(row_k31())
    except ct.Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return EXIT_REFUSED
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("row")
    d = sub.add_parser("dryrun")
    d.add_argument("--out", required=True)
    t = sub.add_parser("start")
    t.add_argument("--out", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", required=True)
    r.add_argument("--seeds", default=None)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--run-dir", default=str(k18.RUN_DIR))
    r.add_argument("--decision", default=str(k18.DECISION))
    r.add_argument("--manifest", default=str(MANIFEST))
    e = sub.add_parser("read")
    e.add_argument("--runs", required=True)
    e.add_argument("--out", default=str(RESULTS))
    a = ap.parse_args(argv)
    return {"row": cmd_row, "dryrun": cmd_dryrun, "start": cmd_start, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
