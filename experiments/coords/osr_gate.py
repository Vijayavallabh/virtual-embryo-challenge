#!/usr/bin/env python
"""The OSR gate: target-free, scoreless. From released input stages only, how large OSR's one geometric step (a
stretch of the overlay union along its own principal axes towards the calendar-fraction log interpolation of the two
brackets' normalised principal SDs) would be, at the test pair and on the nine practice tasks. It computes nothing
against any target and no score, prints no cell count, no board limit and no score.

  pair --name NAME --data DIR --out DIR --commit SHA
      NAME is `test` (E7.5 and E8.0, experiments/coords/e49_shape_interp_read.py's DESCRIBE_PAIR, with the panel and
      count rule of T2:embryo:val_interp, t = 1/2) or a task of experiments/harness/nine_tasks.py (its before and after
      inputs only, with its t, its board panel and count rule; its target is never opened, loaded or stat'ed). Refused
      (exit 2) unless HEAD is SHA and the tracked tree is clean. For each build seed 1,000,001 to 1,000,010 it runs the
      specification's computation and writes OUT/NAME.json; OUT/NAME.exit holds "running", then the exit code: 0, 3 if
      any seed raised or was refused (recorded for that seed: stop clause S5), 1 an error outside the seeds, 2 refused.
  read --out DIR
      reads the ten pair files and applies the stop clauses S1 to S5 into OUT/gate.json and OUT/gate.txt (verdict PASS
      or STOP); exit 0 when both files are written, whatever the verdict.

Specification (the gate's `computation`, word for word): res = interpolate_clouds(Cb, Ca, Eb, Ea, t, n, seed,
target_rms=1.0, t_pos=0.0, positions=POSITIONS[1], shape_form='overlay'); U = res.coords (float64); Zb, Za = the e49
script's _frame of each bracket; blocks = bracket_blocks(Zb, Za); target = log_interp_sds(normalised_sds(Zb),
normalised_sds(Za), t); s_o, Vt = principal(U); f = stretch_factors(s_o, target, blocks); with blocks 1|2,3
g = log(f[0] / f[1]); for any other blocks g = max over axes of |log(f_j / f_k)| with the sign of log(f[0] / f[2]).
Recorded, deciding nothing: g_anti (the same with anti=True), the thin SD of U and of the stretched union against both
brackets', the angle in degrees between Vt[0] and the coordinate z axis of the earlier frame. Per pair, the median g.

Readings where the specification left a choice (declared):
  * count on a twin-board task: the harness's synthetic board takes max_cells = max(5000, 3 max(truth, ref)), and the
    truth term needs the target's row count, which this gate may not read. The target term is left out:
    max_cells = harness.pseudo_counts(0, n_ref)["max_cells"] = max(5000, 3 ref) (n_ref from the reference file as
    released, as run_task reads it). It equals the harness's count whenever the target's truth half is no larger than
    the reference's 10% copy. A validation-board task (e75) takes its published board, as task_board does.
  * the stretched union: mean(U) + (U - mean(U)) @ stretch_matrix(Vt, f); "against both brackets'": the thin SD of
    normalised_sds(Zb) and of normalised_sds(Za) are recorded beside it.
  * the angle: degrees(arccos(|Vt[0] . e_z|)), in [0, 90] (a principal axis has no sign); the earlier frame differs
    from U's coordinates by a shift and a scale only, so its z axis is e_z.
  * g's sign is np.sign(log(f[0] / f[2])), so g = 0 when that log is 0 (the identity blocks 1,2,3 give g = 0).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import traceback
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
import e49_shape_interp_read as e49  # noqa: E402  (it also puts experiments/harness on the path)
import nine_tasks as nt  # noqa: E402

from vec import harness  # noqa: E402
from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402

SCRIPT = "experiments/coords/osr_gate.py"
TEST = "test"
TEST_LIKE = ("half", "late", "new1", "new2")
TASK_NAMES = ("half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3")
NAMES = (TEST, *TASK_NAMES)
SEEDS = e49.DESCRIBE_SEEDS  # 1,000,001 to 1,000,010
TEST_T = 0.5
SHAPE_FORM = "overlay"
TARGET_RMS = 1.0
T_POS = 0.0
SMALL = 0.05  # S1 and S4: |g| below this is "within noise"
WANT_BLOCKS = "1|2,3"
REPORTED_ONLY = ("If the test-pair |median g| lies outside the range of the four test-like tasks' |median g|, any "
                 "practice gain is shown at the test as a range from zero.")


# ------------------------------------------------------------------------------------------------ hooks (tests swap)
def _load(path: Path, genes: list[str]):
    from vec.prep import load_stage

    return load_stage(path, genes=genes, drop_duplicates=True)


def _rows_as_released(path: Path, genes: list[str]) -> int:
    from vec.io import read_stage

    return int(read_stage(path, genes).n_obs)


def _boards():
    from vec.boards import load_boards

    return load_boards()


def _tasks() -> dict:
    return nt.all_tasks()


# ------------------------------------------------------------------------------------------------ checks
def commit_problems(commit: str, repo: Path = REPO) -> list[str]:
    """Why this checkout may not run the gate at `commit` (empty: it may): HEAD and a clean tracked tree."""
    env = {**os.environ, "GIT_OPTIONAL_LOCKS": "0"}
    head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, env=env)
    st = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"],
                        capture_output=True, text=True, env=env)
    why = []
    if head.returncode != 0 or head.stdout.strip() != commit:
        why.append(f"HEAD is {head.stdout.strip()[:12] or 'unknown'}, not {commit[:12]}")
    if st.returncode != 0 or st.stdout.strip():
        why.append("the tracked tree is not clean")
    return why


def own_code_problems() -> list[str]:
    src = (REPO / "src").resolve()
    return [f"{m.__name__} is loaded from {m.__file__}, not this checkout's src/" for m in (harness, si, tb)
            if src not in Path(m.__file__).resolve().parents]


# ------------------------------------------------------------------------------------------------ the computation
def g_of(f, blocks) -> float:
    """The specification's g from the stretch factors f and the blocks."""
    f = np.asarray(f, np.float64)
    if si.blocks_text(blocks) == WANT_BLOCKS:
        return float(np.log(f[0] / f[1]))
    big = max(abs(float(np.log(f[j] / f[k]))) for j in range(3) for k in range(3))
    return float(np.sign(np.log(f[0] / f[2])) * big)


def seed_quantities(Cb, Ca, Eb, Ea, t: float, n: int, seed: int) -> dict:
    """One pair-seed: the specification's computation, word for word, and the recorded-only quantities. No count."""
    res = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, n, seed, target_rms=TARGET_RMS, t_pos=T_POS,
                                positions=tb.POSITIONS[1], shape_form=SHAPE_FORM)
    U = np.asarray(res.coords, dtype=np.float64)
    Zb, Za = e49._frame(Cb), e49._frame(Ca)
    blocks = si.bracket_blocks(Zb, Za)
    s_b, s_a = si.normalised_sds(Zb), si.normalised_sds(Za)
    target = si.log_interp_sds(s_b, s_a, t)
    s_o, Vt = si.principal(U)
    f = si.stretch_factors(s_o, target, blocks)
    g = g_of(f, blocks)
    target_anti = si.log_interp_sds(s_b, s_a, t, anti=True)
    f_anti = si.stretch_factors(s_o, target_anti, blocks)
    m = U.mean(0)
    stretched = m + (U - m) @ si.stretch_matrix(Vt, f)
    angle = float(np.degrees(np.arccos(min(1.0, abs(float(Vt[0][2]))))))
    return {"blocks": si.blocks_text(blocks), "g": g, "g_anti": g_of(f_anti, blocks),
            "f": [float(x) for x in f], "f_anti": [float(x) for x in f_anti],
            "target_sds": [float(x) for x in target], "target_sds_anti": [float(x) for x in target_anti],
            "own_sds": [float(x) for x in s_o],
            "thin_sd_union": float(s_o[-1]), "thin_sd_stretched": float(si.normalised_sds(stretched)[-1]),
            "thin_sd_before": float(s_b[-1]), "thin_sd_after": float(s_a[-1]),
            "angle_axis1_z_deg": angle}


def pair_inputs(name: str, data: Path) -> dict:
    """The pair's brackets (coordinates, panel expression), t, count and what was read. Never a target."""
    from vec.io import COORD_KEY, dense

    boards = _boards()
    if name == TEST:
        board_key, files, t = e49.DESCRIBE_BOARD, dict(e49.DESCRIBE_PAIR), TEST_T
        board = boards[board_key]
        count_rule = f"harness.count_for({board_key}, before's rows after drop_duplicates)"
    else:
        task = _tasks()[name]
        board_key, files = task.panel_board, {r: task.inputs[r] for r in ("before", "after")}
        view = harness.TaskView.of(task)
        lo, hi = view.role_days["before"], view.role_days["after"]
        t = (view.target_day - lo) / (hi - lo)  # as harness._fgw_* reads it: STAGE_DAYS, no file
        if task.validation_board:
            board = replace(boards[task.validation_board], anchors={})
            count_rule = f"harness.count_for({task.validation_board} as published, before's rows)"
        else:
            twin = boards[task.twin_board]
            n_ref = _rows_as_released(data / task.reference, list(twin.genes))
            pc = harness.pseudo_counts(0, n_ref)  # the target term left out (module docstring)
            board = replace(twin, key=f"pseudo:{task.key}", split="pseudo", anchors={}, max_cells=pc["max_cells"],
                            min_cells=1, truth_cells=None, ref_cells=pc["ref_cells"])
            count_rule = (f"harness.count_for(pseudo board on {task.twin_board}, before's rows), max_cells = "
                          "max(5000, 3 ref_cells): the target term left out")
    if e49._is_test_board(board) or board.counts_hidden:
        raise PermissionError("the gate reads a validation or pseudo board's panel and count rule only")
    genes = list(board.genes)
    st = {r: _load(data / f, genes) for r, f in files.items()}
    for r, a in st.items():
        missing = [g for g in genes if g not in a.var_names]
        if missing:
            raise ValueError(f"{name}: input {r} lacks {len(missing)} panel genes")
    Cb, Ca = (np.asarray(st[r].obsm[COORD_KEY], float)[:, :3] for r in ("before", "after"))
    Eb, Ea = (dense(st[r][:, genes].X).astype(np.float32) for r in ("before", "after"))
    n = harness.count_for(board, st["before"].n_obs)
    return {"Cb": Cb, "Ca": Ca, "Eb": Eb, "Ea": Ea, "t": float(t), "n": n, "files": files,
            "panel_board": board_key, "count_rule": count_rule}


def run_pair(name: str, data: Path, out: Path, commit: str) -> int:
    """Write out/<name>.json; return 0, or 3 if any seed raised or was refused, or 1 on an error outside the seeds."""
    rec: dict = {"name": name, "commit": commit, "script": SCRIPT, "seeds": list(SEEDS), "shape_form": SHAPE_FORM}
    code = 0
    try:
        p = pair_inputs(name, Path(data))
        rec.update(t=p["t"], files=p["files"], panel_board=p["panel_board"], count_rule=p["count_rule"])
        per_seed = {}
        for seed in SEEDS:
            try:
                per_seed[str(seed)] = seed_quantities(p["Cb"], p["Ca"], p["Eb"], p["Ea"], p["t"], p["n"], seed)
            except Exception as e:  # stop clause S5: recorded for the seed, the run goes on
                per_seed[str(seed)] = {"error": f"{type(e).__name__}: {e}"}
                code = 3
        rec["per_seed"] = per_seed
        gs = [v["g"] for v in per_seed.values() if "g" in v]
        rec["median_g"] = float(np.median(gs)) if gs else None
        rec["median_g_over"] = len(gs)
    except Exception as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["traceback"] = traceback.format_exc(limit=5)
        code = 1
    (out / f"{name}.json").write_text(json.dumps(rec, indent=1))
    return code


def cmd_pair(a) -> int:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / f"{a.name}.json").exists():
        print(f"refused: {out / (a.name + '.json')} exists", file=sys.stderr)
        return 2
    ex = out / f"{a.name}.exit"
    ex.write_text("running\n")
    why = commit_problems(a.commit) + own_code_problems()
    if why:
        print("refused: " + "; ".join(why), file=sys.stderr)
        code = 2
    else:
        code = run_pair(a.name, Path(a.data), out, a.commit)
    ex.write_text(f"{code}\n")
    return code


# ------------------------------------------------------------------------------------------------ the reading
def _sign(x: float) -> int:
    return int(np.sign(x))


def reading(out: Path) -> dict:
    """The specification's stop_result applied to the ten pair files in out."""
    s5, pairs, commits = [], {}, set()
    for name in NAMES:
        p = Path(out) / f"{name}.json"
        try:
            r = json.loads(p.read_text())
        except Exception as e:
            s5.append(f"{name}: no readable pair file ({type(e).__name__})")
            continue
        commits.add(r.get("commit"))
        if "error" in r:
            s5.append(f"{name}: {r['error']}")
            continue
        per = r.get("per_seed") or {}
        if sorted(per) != sorted(str(s) for s in SEEDS):
            s5.append(f"{name}: seeds are not {SEEDS[0]} to {SEEDS[-1]}")
        for s, v in per.items():
            if "error" in v or "g" not in v or v["g"] is None or not math.isfinite(v["g"]):
                s5.append(f"{name} seed {s}: {v.get('error', 'no finite g')}")
        pairs[name] = r
    med = {n: pairs[n]["median_g"] for n in NAMES if n in pairs and pairs[n].get("median_g") is not None}
    held: dict[str, str] = {}
    notes: dict[str, str] = {}
    tp = pairs.get(TEST)
    if tp is not None:
        gs = [v["g"] for v in tp["per_seed"].values() if "g" in v and v["g"] is not None]
        small = sum(abs(g) < SMALL for g in gs)
        if small >= 8:
            held["S1"] = f"|g| < {SMALL} in {small} of 10 test-pair seeds"
        signs = [_sign(g) for g in gs] + [None] * (len(SEEDS) - len(gs))  # a seed without g agrees with no sign
        top = max((signs.count(s) for s in (-1, 0, 1)), default=0)
        if len(signs) - top > 2:
            held["S2"] = f"sign(g) differs from the majority sign in {len(signs) - top} of 10 test-pair seeds"
        blocks = sorted({v["blocks"] for v in tp["per_seed"].values() if "blocks" in v})
        if blocks != [WANT_BLOCKS]:
            held["S3"] = f"test-pair blocks are {blocks or 'unknown'}, not {WANT_BLOCKS}"
    else:
        notes["S1-S3"] = "not evaluable: no test-pair file"
    like = {n: med[n] for n in TEST_LIKE if n in med}
    if TEST in med and len(like) == len(TEST_LIKE):
        opp = [n for n, m in like.items() if _sign(med[TEST]) * _sign(m) < 0]
        if len(opp) >= 3:
            held["S4"] = f"test-pair median g has the opposite sign to the median g of {', '.join(opp)}"
        elif all(abs(m) < SMALL for m in like.values()):
            held["S4"] = f"all four test-like tasks have |median g| < {SMALL}"
    else:
        notes["S4"] = "not evaluable: a test-pair or test-like median g is missing"
    if s5:
        held["S5"] = "; ".join(s5)
    inside = None
    if TEST in med and len(like) == len(TEST_LIKE):
        lo, hi = min(abs(m) for m in like.values()), max(abs(m) for m in like.values())
        inside = bool(lo <= abs(med[TEST]) <= hi)
    clauses = [c for c in ("S1", "S2", "S3", "S4", "S5") if c in held]
    return {"verdict": "STOP" if clauses else "PASS", "stopped_by": clauses, "why": held, "not_evaluable": notes,
            "median_g": med, "commits": sorted(str(c) for c in commits),
            "reported_only": {"test_abs_median_g_inside_test_like_range": inside, "statement": REPORTED_ONLY},
            "script": SCRIPT}


def text(res: dict) -> str:
    lines = [f"OSR gate: {res['verdict']}" + (f" ({', '.join(res['stopped_by'])})" if res["stopped_by"] else "")]
    for c in res["stopped_by"]:
        lines.append(f"  {c}: {res['why'][c]}")
    for c, v in res["not_evaluable"].items():
        lines.append(f"  {c}: {v}")
    lines.append("median g per pair:")
    for n in NAMES:
        m = res["median_g"].get(n)
        lines.append(f"  {n:<5} " + ("missing" if m is None else f"{m:+.4f}"))
    ins = res["reported_only"]["test_abs_median_g_inside_test_like_range"]
    lines.append("reported only, deciding nothing: test-pair |median g| inside the test-like range: "
                 + ("not evaluable" if ins is None else str(ins)))
    lines.append("  " + res["reported_only"]["statement"])
    return "\n".join(lines) + "\n"


def cmd_read(a) -> int:
    out = Path(a.out)
    res = reading(out)
    (out / "gate.json").write_text(json.dumps(res, indent=1))
    (out / "gate.txt").write_text(text(res))
    print(text(res), end="")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pair")
    p.add_argument("--name", choices=NAMES, required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True)
    r = sub.add_parser("read")
    r.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    return {"pair": cmd_pair, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
