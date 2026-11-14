#!/usr/bin/env python
"""Row N51 of experiments/harness/TRIES.md (a diagnostic, not a try; target-free): is the entry's coupling of the
misshapen kind on each draw of the nine practice tasks?

The embryo entry `fgw` moves cells along an FGW coupling of the two bracket stages. On some draws the coupling's
partners sit mirrored across one direction and the moved cloud comes out flattened. The diagnostic of 2026-11-03
(`diag_emb_coupling.py`, row E10's gate) measured this for the two counted backtests and the test pair. This script
measures it on all nine practice tasks of `nine_tasks.py` (`all_tasks()`), at the harness's ten build seeds
(`vec.harness.BUILD_SEED_OFFSET + s`, s in 0..9), with the diagnostic's own two measures, imported and not copied:
  xcorr_min    the smallest eigenvalue of the symmetrised whitened cross-correlation of the anchors and their
               partners (`diag_emb_coupling._measure`);
  thin_ratio   the moved cloud's thinnest principal SD over its RMS radius, against the thinner bracket's;
and the diagnostic's committed flag, defect = xcorr_min < DEFECT_XCORR or thin_ratio < DEFECT_THIN (0.3 and 0.9, no
new number), with `by`: which of the two conditions set it ("xcorr", "thin", "both"; null when not flagged).

The coupling is the harness's (`vec.harness._fgw`): the task's two bracket files (its "before" and "after" inputs)
loaded as `vec.harness.run_task` loads them (`vec.prep.load_stage` on the task's panel board's genes, duplicates
dropped), the panel of the task's panel board (`PseudoTask.panel_board`), the fraction the harness computes
((target day - before day) / (after day - before day), from `vec.harness.TaskView`), the entry's moves ("fgw"), the
build seed. The draw's size and the size rule change neither the coupling nor the moves the measures read, so the
draw uses the diagnostic's N_DRAW and target_rms = 1 (the harness's count rule reads the target's size, which this
script never touches).

Target-free: it opens the two brackets of each task and nothing else; a path equal to the task's target is refused
before it is opened. It scores nothing and writes no cell count of any stage.

Reproduction gate, checked before n51.json and n51.txt are written and, in a run, before any part file is written (the
counted backtests among the run's tasks are measured first): on the two counted backtests (e75, e725), xcorr_min and
thin_ratio equal the committed ones in results/diag_emb_coupling/diag_emb_coupling.json within GATE_TOL on every seed;
otherwise the script stops and writes nothing.

Refused: a run without row N51 in TRIES.md, a dirty tree, an existing output.

  run   [--tasks NAME ...] --out DIR   writes n51_<task>.json per task; without --tasks all nine, then the merge
  merge --out DIR                      applies the gate to the nine part files and writes n51.json and n51.txt

Where it runs. Some builds differ between the two machines (docs/COMPUTE.md, "What reproduces where"), and the flag
must come from the coupling on the machine that produced the committed arm scores of rows N45-N47: run it on the DGX
campaign checkout (~/v/vec-campaign, where the campaign's harness runs go; AGENTS.md rule 8), with the BLAS threads
capped at 16 as below. The gate ties this run to the diagnostic's committed run only, not to the score runs.

Usage (DGX campaign checkout, repo root):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
    uv run --frozen --no-sync python -I -B experiments/harness/n51_coupling_defect_nine.py run \
    --out experiments/harness/results/n51_coupling_defect_nine
"""

from __future__ import annotations

import argparse
import json
import sys
import textwrap
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))
import diag_emb_coupling as dc  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import harness, t2_build  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.prep import load_stage  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW = "N51"
TRIES = HERE / "TRIES.md"
OUT = HERE / "results" / "n51_coupling_defect_nine"
COMMITTED = HERE / "results" / "diag_emb_coupling" / "diag_emb_coupling.json"
SEEDS = dc.SEEDS  # s in 0..9, build seed harness.BUILD_SEED_OFFSET + s, as the harness's
GATE_TASKS = ("e75", "e725")  # the two counted embryo backtests; the diagnostic's pairs under the same task keys
GATE_TOL = 1e-6
MEASURES = ("xcorr_min", "thin_ratio")


def nine() -> dict:
    """The nine practice tasks, by name (nine_tasks.all_tasks)."""
    return nt.all_tasks()


def _genes(t) -> list[str]:
    """The task's panel, as vec.harness.run_task takes it."""
    return list(load_boards()[t.panel_board].genes)


def _read(path: Path, genes: list[str]):
    """A bracket stage, as vec.harness.run_task loads a T2 input."""
    return load_stage(path, genes=genes, drop_duplicates=True)


def _state() -> tuple[str, bool, str]:
    s = repo_state(REPO)
    return s.commit, s.unclean == (), s.utc


def bracket_paths(t, data_dir: Path) -> dict[str, Path]:
    """The task's two brackets; refuses a path that is the task's target."""
    paths = {r: Path(data_dir) / t.inputs[r] for r in ("before", "after")}
    target = Path(data_dir) / t.target
    if any(p == target for p in paths.values()):
        raise harness.LeakageError(f"{t.key}: a bracket path is the task's target")
    return paths


def fraction(t) -> float:
    """The target's time fraction of the bracket, as vec.harness._fgw computes it."""
    v = harness.TaskView.of(t)
    lo, hi = v.role_days["before"], v.role_days["after"]
    return (v.target_day - lo) / (hi - lo)


def classify(xcorr_min: float, thin_ratio: float) -> tuple[bool, str | None]:
    """The diagnostic's committed flag and which condition set it."""
    x, th = bool(xcorr_min < dc.DEFECT_XCORR), bool(thin_ratio < dc.DEFECT_THIN)
    return x or th, ("both" if x and th else "xcorr" if x else "thin" if th else None)


def measure_draw(inputs: tuple, frac: float, seed: int) -> dict[str, float]:
    """The entry's coupling at this build seed, and the diagnostic's two measures for its moves."""
    Cb, Ca, Eb, Ea, lab_b, lab_a = inputs
    res = t2_build.interpolate_clouds(Cb, Ca, Eb, Ea, frac, dc.N_DRAW, seed, target_rms=1.0, with_coupling=True)
    m = dc._measure(res, frac, lab_b, lab_a, dc._nn_rows(seed, len(Cb), len(Ca)))
    return {k: float(m[k]) for k in MEASURES}


def measure_task(t, data_dir: Path, cache: dict | None = None, log=print) -> dict:
    """One task's record: brackets, panel board, fraction and the ten draws."""
    cache = {} if cache is None else cache
    genes = _genes(t)
    paths = bracket_paths(t, data_dir)
    stages = []
    for r in ("before", "after"):
        k = (str(paths[r]), tuple(genes))
        if k not in cache:
            cache[k] = _read(paths[r], genes)
        stages.append(cache[k])
    inputs = dc._inputs(*stages, genes)
    frac = fraction(t)
    draws = []
    for s in SEEDS:
        t0 = time.time()
        m = measure_draw(inputs, frac, harness.BUILD_SEED_OFFSET + s)
        defect, by = classify(m["xcorr_min"], m["thin_ratio"])
        draws.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, **m, "defect": defect, "by": by})
        log(f"{t.key} seed {s}: xcorr_min {m['xcorr_min']:.3f} thin_ratio {m['thin_ratio']:.3f} defect {defect}"
            f"{f' ({by})' if by else ''} ({time.time() - t0:.1f} s)", flush=True)
    return {"key": t.key, "brackets": [harness.stem(t.inputs["before"]), harness.stem(t.inputs["after"])],
            "panel_board": t.panel_board, "fraction": frac, "draws": draws,
            "n_defect": sum(d["defect"] for d in draws)}


def gate(records: dict, committed: dict) -> tuple[float | None, list[str]]:
    """The largest absolute difference between the records' measures and the committed diagnostic's on the counted
    backtests among `records` (task name -> record), and the reasons the gate fails."""
    worst, why = None, []
    for name in GATE_TASKS:
        if name not in records:
            continue
        rec = records[name]
        pair = committed.get("pairs", {}).get(rec["key"])
        if pair is None:
            why.append(f"{name}: the committed diagnostic has no pair {rec['key']}")
            continue
        if (pair["before"], pair["after"]) != tuple(rec["brackets"]):
            why.append(f"{name}: the committed pair's brackets are {pair['before']}, {pair['after']}")
        old = {d["seed"]: d["fgw"] for d in pair["draws"]}
        for d in rec["draws"]:
            if d["seed"] not in old:
                why.append(f"{name}: the committed diagnostic has no seed {d['seed']}")
                continue
            for k in MEASURES:
                diff = abs(float(d[k]) - float(old[d["seed"]][k]))
                worst = diff if worst is None else max(worst, diff)
                if not diff <= GATE_TOL:
                    why.append(f"{name} seed {d['seed']}: {k} differs from the committed one by {diff:.2e}")
    return worst, why


def _committed() -> dict:
    return json.loads(COMMITTED.read_text())


def _guard(out_files) -> tuple[str, bool, str]:
    if ROW not in n12.tries_rows(TRIES):
        sys.exit(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    existing = [str(p) for p in out_files if Path(p).exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    commit, clean, utc = _state()
    if not clean:
        sys.exit(f"the tree is not clean: {ROW} runs on a committed tree")
    return commit, clean, utc


def part_path(out: Path, name: str) -> Path:
    return Path(out) / f"n51_{name}.json"


def cmd_run(a) -> None:
    tasks = nine()
    names = list(a.tasks) if a.tasks else list(tasks)
    unknown = [n for n in names if n not in tasks]
    if unknown:
        sys.exit(f"unknown tasks {unknown}; known: {sorted(tasks)}")
    out = Path(a.out)
    files = [part_path(out, n) for n in names]
    if not a.tasks:
        files += [out / "n51.json", out / "n51.txt"]
    commit, clean, _ = _guard(files)
    committed = _committed() if any(n in GATE_TASKS for n in names) else None
    cache: dict = {}
    # the counted backtests are measured and gated first, so that a failed gate leaves no part file behind
    gated = {}
    for n in [n for n in names if n in GATE_TASKS]:
        rec = measure_task(tasks[n], Path(a.data), cache)
        worst, why = gate({n: rec}, committed)
        if why:
            sys.exit("reproduction gate failed, nothing written: " + "; ".join(why))
        rec["gate_max_abs_difference"] = worst
        gated[n] = rec
    out.mkdir(parents=True, exist_ok=True)
    for n in [*gated, *(n for n in names if n not in GATE_TASKS)]:
        rec = gated[n] if n in gated else measure_task(tasks[n], Path(a.data), cache)
        _, clean_end, utc = _state()
        rec.update({"task": n, "row": ROW, "commit": commit, "clean": bool(clean and clean_end), "utc": utc})
        part_path(out, n).write_text(json.dumps(rec, indent=1) + "\n")
        print(f"wrote {part_path(out, n)}", flush=True)
    if not a.tasks:
        merge(out)


def merge(out: Path) -> dict:
    """Reads the nine part files, applies the gate and writes n51.json and n51.txt."""
    out = Path(out)
    tasks = nine()
    commit, clean, utc = _guard([out / "n51.json", out / "n51.txt"])
    missing = [n for n in tasks if not part_path(out, n).exists()]
    if missing:
        sys.exit(f"no part file for {missing}")
    parts = {n: json.loads(part_path(out, n).read_text()) for n in tasks}
    bad = [n for n, p in parts.items() if p["key"] != tasks[n].key or [d["seed"] for d in p["draws"]] != list(SEEDS)]
    if bad:
        sys.exit(f"part files with another task key or other seeds: {bad}")
    if {p["commit"] for p in parts.values()} != {commit}:
        sys.exit(f"the part files' commits are not this tree's ({commit})")
    if not all(p["clean"] for p in parts.values()):
        sys.exit("a part file comes from a dirty tree")
    worst, why = gate(parts, _committed())
    if why or worst is None:
        sys.exit("reproduction gate failed, nothing written: " + ("; ".join(why) or "no counted backtest measured"))
    keep = ("key", "brackets", "panel_board", "fraction", "n_defect", "draws")
    report = {"what": __doc__.split("\n")[0], "row": ROW, "commit": commit, "clean": clean, "utc": utc,
              "thresholds": {"DEFECT_XCORR": dc.DEFECT_XCORR, "DEFECT_THIN": dc.DEFECT_THIN},
              "seeds": list(SEEDS), "build_seed_offset": harness.BUILD_SEED_OFFSET,
              "gate": {"tasks": list(GATE_TASKS), "tol": GATE_TOL, "max_abs_difference": worst,
                       "committed": str(COMMITTED.relative_to(REPO)), "committed_commit": _committed().get("commit")},
              "parts_utc": {n: p["utc"] for n, p in parts.items()},
              "tasks": {n: {k: p[k] for k in keep} for n, p in parts.items()}}
    txt = text(report)
    (out / "n51.json").write_text(json.dumps(report, indent=1) + "\n")
    (out / "n51.txt").write_text(txt)
    print(txt, end="")
    return report


def text(report: dict) -> str:
    lines = textwrap.wrap(
        f"Row {report['row']} (a diagnostic, not a try; target-free): the entry's FGW coupling on the nine practice "
        f"tasks at build seeds {report['build_seed_offset']:,} + 0..9, flagged by the committed thresholds of "
        f"diag_emb_coupling.py (xcorr_min < {report['thresholds']['DEFECT_XCORR']} or thin_ratio < "
        f"{report['thresholds']['DEFECT_THIN']}). Commit {report['commit']}, clean tree: {report['clean']}, "
        f"{report['utc']}. Reproduction gate on {', '.join(report['gate']['tasks'])}: largest difference "
        f"{report['gate']['max_abs_difference']:.1e} (tolerance {report['gate']['tol']:.0e}).", 118)
    lines += ["", "  task  brackets              t      flagged  per seed 0..9 (x: xcorr, t: thin, b: both, .: none)"]
    mark = {"xcorr": "x", "thin": "t", "both": "b", None: "."}
    for n, r in report["tasks"].items():
        lines.append(f"  {n:<5} {' -> '.join(r['brackets']):<21} {r['fraction']:.3f}  {r['n_defect']:>2}/10   "
                     + " ".join(mark[d["by"]] for d in r["draws"]))
    lines += ["", "  xcorr_min per seed / thin_ratio per seed:"]
    for n, r in report["tasks"].items():
        lines.append(f"  {n:<5} " + " ".join(f"{d['xcorr_min']:+.2f}" for d in r["draws"]))
        lines.append(f"  {'':<5} " + " ".join(f"{d['thin_ratio']:5.2f}" for d in r["draws"]))
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--tasks", nargs="+", default=None, help="a subset of the nine (one part file each); no merge")
    r.add_argument("--out", default=str(OUT))
    r.add_argument("--data", default=str(harness.DATA_DIR))
    m = sub.add_parser("merge")
    m.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    if a.cmd == "run":
        cmd_run(a)
    else:
        merge(Path(a.out))


if __name__ == "__main__":
    main()
