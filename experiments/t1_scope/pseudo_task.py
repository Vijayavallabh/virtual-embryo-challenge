"""T1 scope pseudo-task: chamber re-weighting of real reference-stage cells, and a damped per-type count-space shift
on top of it, against copy_last (docs/PLAN.md section 4.1; screened digest gap-fill 1, which ran the same geometry).

Geometry (released T1 stages only): reference and copy_last = E8.5; truth = one E9.5 library cut to the chamber
categories of configs/t1_chamber_classes.yaml plus a small surrounding contamination of foregut-type cells (5% of the
truth by default, standing in for the few oesophagus-junction cells of the E12.5 dissection). Scored with the T1 board
emulator (vec.evaluate: 10% working copies, truth half against the other half as ceiling, the reference working copy
as floor), anchors computed per seed; every candidate is rebuilt per seed from its own sampling stream.

Candidates (n_cells real E8.5 cells each; the reweighting candidates go through vec.t1_scope.reweight):
  copy_last         a uniform sample of E8.5
  oracle            re-weighted to the truth's own category shares; surrounding drawn from the contamination labels
  oracle_any_surr   the same with surrounding drawn from every E8.5 surrounding label
  cm_minus/cm_plus  oracle with the CM share of heart cells moved by -/+ cm_delta (10 points by default)
  drop_surr_only    E8.5's own composition without the surrounding category (adjacent types kept)
  blood_plus        oracle with blood at +10 points of all cells (other categories scaled)
  shift_<d>         oracle + momentum_shift(start=E8.5, ref=the other E9.5 library, damping d). No T1 stage precedes
                    E8.5, so honest momentum does not exist in this geometry: the change is measured on the truth
                    stage's other library, an oracle direction (upper bound for a per-type count-space shift)
  flip_<d>          the same change with the sign reversed (the cost of a wrong direction)
  momentum_<d>      only with --momentum-start: oracle + momentum_shift(start=that stage, ref=the reference stage,
                    damping d), i.e. real momentum from the interval before the reference
  oracle_no_epi     ablation: oracle with the epicardium share redistributed (no JCF stand-ins at E8.5)
  oracle_n<k>       ablation: oracle with k cells (--small-n, 2000 by default)
--only limits the run to the named candidates (copy_last and oracle are always kept for the paired differences).
E8.5 has no endocardium, vascular or epicardium label of its own, so the config's E8.5 rules apply. Its override
(Endothelium as endocardium) is a remap: those cells count as endocardium's own, in the draws and in the momentum.
Its stand-ins (Endothelium as vascular endothelium, JCF as epicardium) fill re-weighting quotas; a type drawn from
stand-ins at either stage gets no momentum (reported as skipped) unless --standin-momentum is given
(standin_momentum=True in vec.t1_scope.momentum_shift).

After the 2026-10-20 release the same script runs the geometry where momentum is honest (map the new labels in
configs/t1_chamber_classes.yaml first; classify refuses unmapped labels):
  --ref data/raw/T1/E9.5_RNA.h5ad --truth-stage data/raw/T1/E10.5_RNA.h5ad --momentum-start data/raw/T1/E8.5_RNA.h5ad
The truth library needs cells of the --contamination-labels, and each momentum endpoint (sibling libraries, reference,
momentum start) at least --min-cells of them in all, as momentum_shift needs for surrounding; the counts are checked
at load time and the run refuses otherwise.

Usage (DGX, repo root; ~10 min with 5 jobs; run in tmux with a log under runs/logs/):
  uv run python experiments/t1_scope/pseudo_task.py --ref data/raw/T1/E8.5_RNA.h5ad \
      --truth-stage data/raw/T1/E9.5_RNA.h5ad --truth-library 2 --seeds 0,1,2,3,4,5,6,7,8,9 --jobs 5
Writes aggregates only (scores per seed and candidate, summaries, compositions, shift factors per type, run record)
to experiments/t1_scope/results/<tag>/, all at the end and only if the procedure files still hash as they did at the
start (run.json last). Nothing per cell leaves runs/.
"""

from __future__ import annotations

import hashlib
import json
import multiprocessing as mp
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import typer
from rich.console import Console
from rich.table import Table

from vec.boards import RANKED_METRICS, REPO_ROOT, get_board
from vec.evaluate import (
    DEFAULT_BUILD,
    SUBSAMPLE_FRAC,
    Stage,
    aggregate,
    panel_metrics,
    score_board,
    skills_from_anchors,
    split_half,
    subsample,
)
from vec.prep import load_stage
from vec.provenance import stamp_path
from vec.t1_scope import (
    build_submission,
    classify,
    composition,
    load_shares,
    load_spec,
    momentum_shift,
    reweight,
    share_of,
    with_share,
    write_checked,
)

app = typer.Typer(add_completion=False)
console = Console(emoji=False, highlight=False)
RESULTS = Path(__file__).resolve().parent / "results"
TASK = "T1"
EXTRA_METRICS = ("energy_distance", "pb_rel_err", "library_size_ratio", "variance_ratio", "pseudobulk_pearson",
                 "_de_raw", "_de_chance", "_n_up", "_n_dn")
_DATA: dict = {}  # set in the parent before the worker pool forks


def sha256_file(path: Path, chunk: int = 1 << 23) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


PROCEDURE_FILES = ("src/vec/t1_scope.py", "configs/t1_chamber_classes.yaml", "experiments/t1_scope/pseudo_task.py",
                   "experiments/t1_direction/results/held_flat_genes.csv")


def git_state() -> dict:
    """HEAD, the procedure files that differ from it (untracked files listed one by one) and their sha256."""

    def run(*args):
        try:
            return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unknown"

    changed = {ln[3:] for ln in run("status", "--porcelain", "--untracked-files=all").splitlines()}
    return {"commit": run("rev-parse", "--short", "HEAD"),
            "uncommitted": [f for f in PROCEDURE_FILES if f in changed],
            "sha256": {f: sha256_file(REPO_ROOT / f) for f in PROCEDURE_FILES}}


def write_results(out_dir: Path, files: dict[str, str], git_at_start: dict) -> None:
    """Write the run's files (name -> text) into out_dir, but only if the procedure files still hash as they did at
    the start; otherwise nothing is written. An older run.json is removed first and the new one is written last, so a
    run.json only ever sits beside its own run's outputs."""
    if git_state()["sha256"] != git_at_start["sha256"]:
        raise SystemExit("procedure files changed during the run: nothing was written, since the record would not "
                         "describe the code that ran; rerun")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "run.json").unlink(missing_ok=True)
    for name, text in sorted(files.items(), key=lambda kv: kv[0] == "run.json"):
        (out_dir / name).write_text(text)


def library_of(obs_names) -> np.ndarray:
    """T1 library tag: the obs-name suffix after the last '_' (experiments/t1_direction/t1_stats.py)."""
    return pd.Series(np.asarray(obs_names, dtype=str)).str.split("_").str[-1].to_numpy()


def check_contamination(truth: dict, endpoints: dict, labels: list[str], min_cells: int) -> dict[str, dict[str, int]]:
    """Cells of each contamination label at each stage used (name -> AnnData): the truth library, then the momentum
    endpoints (sibling libraries, reference, momentum start). `within` restricts surrounding to these labels at every
    stage. The truth's contamination is drawn from them, so it needs at least one cell of them; momentum_shift needs
    min_cells of them in all at each endpoint (at least one), or surrounding silently gets no momentum. The run is
    refused at load time when a stage falls short."""
    cells = {name: {lab: int((x.obs["celltype"].astype(str).to_numpy() == lab).sum()) for lab in labels}
             for name, x in {**truth, **endpoints}.items()}
    need = {**dict.fromkeys(truth, 1), **dict.fromkeys(endpoints, max(1, min_cells))}
    short = [f"{name} ({sum(cells[name].values())} < {n})" for name, n in need.items()
             if sum(cells[name].values()) < n]
    if short:
        raise typer.BadParameter(f"too few cells of the contamination labels {labels} at {', '.join(short)}; cells per "
                                 f"stage: {cells}. The truth library needs at least one, and each momentum endpoint "
                                 f"(sibling libraries, reference, momentum start) --min-cells ({min_cells}) in all, "
                                 "or surrounding gets no momentum")
    return cells


def make_truth(stage, library: str, chamber: list[str], contamination: float, contam_labels: list[str], seed: int):
    """One library cut to the chamber categories plus `contamination` (share of the truth) foregut-type cells."""
    lib = library_of(stage.obs_names) == library
    cats = stage.obs["category"].astype(str).to_numpy()
    labels = stage.obs["celltype"].astype(str).to_numpy()
    keep = np.flatnonzero(lib & np.isin(cats, chamber))
    pool = np.flatnonzero(lib & np.isin(labels, contam_labels))
    n_c = int(round(contamination / (1.0 - contamination) * len(keep)))
    if n_c > len(pool):
        raise typer.BadParameter(f"{n_c} contamination cells needed, only {len(pool)} in library {library}")
    pick = np.random.default_rng(seed).choice(pool, n_c, replace=False)
    idx = np.sort(np.concatenate([keep, pick]))
    return stage[idx].copy(), {"library": library, "n_chamber": int(len(keep)), "n_contamination": n_c,
                               "contamination_labels_available": int(len(pool))}


def to_stage(a) -> Stage:
    return Stage(X=a.X, celltype=a.obs["celltype"].astype(str).to_numpy())


def candidates(seed: int) -> tuple[dict, dict, dict]:
    """Every candidate for one seed: name -> Stage, plus re-weighting and shift reports (seed-specific)."""
    d = _DATA
    ref, sib, cfg, mstart = d["ref"], d["sibling"], d["cfg"], d["momentum_start"]
    out, rw_rep, sh_rep = {}, {}, {}
    want = set(cfg["only"]) | {"copy_last", "oracle"} if cfg["only"] else None

    def wanted(name):
        return want is None or name in want

    rng = np.random.default_rng([seed, 99])
    pick = np.sort(rng.choice(ref.n_obs, cfg["n_cells"], replace=False))
    out["copy_last"] = Stage(X=ref.X[pick], celltype=ref.obs["celltype"].astype(str).to_numpy()[pick])
    oracle_pop = None
    for name, (shares, within) in d["share_sets"].items():
        if not wanted(name):
            continue
        pop, rep = reweight(ref, shares, cfg["n_cells"], seed, within=within, fallback=d["fallback"])
        out[name] = to_stage(pop)
        rw_rep[name] = rep.as_dict()
        if name == "oracle":
            oracle_pop = pop
    small = f"oracle_n{cfg['small_n']}"
    if wanted(small):
        shares, within = d["share_sets"]["oracle"]
        pop, rep = reweight(ref, shares, cfg["small_n"], seed, within=within, fallback=d["fallback"])
        out[small] = to_stage(pop)
        rw_rep[small] = rep.as_dict()
    for damping in cfg["dampings"]:
        for sign, prefix in ((1.0, "shift"), (-1.0, "flip")):
            if sign < 0 and damping not in cfg["flip_dampings"] or not wanted(f"{prefix}_{damping:g}"):
                continue
            pop, rep = momentum_shift(ref, sib, oracle_pop, sign * damping, d["held_flat"], within=d["within"],
                                      start_fallback=d["fallback"], ref_fallback=d["sib_fallback"], seed=[seed, 1],
                                      min_cells=cfg["min_cells"], eps=cfg["eps"],
                                      standin_momentum=cfg["standin_momentum"])
            name = f"{prefix}_{damping:g}"
            out[name] = to_stage(pop)
            sh_rep[name] = rep.as_dict()
        if mstart is not None and wanted(f"momentum_{damping:g}"):
            pop, rep = momentum_shift(mstart, ref, oracle_pop, damping, d["held_flat"], within=d["within"],
                                      start_fallback=d["mstart_fallback"], ref_fallback=d["fallback"], seed=[seed, 2],
                                      min_cells=cfg["min_cells"], eps=cfg["eps"],
                                      standin_momentum=cfg["standin_momentum"])
            name = f"momentum_{damping:g}"
            out[name] = to_stage(pop)
            sh_rep[name] = rep.as_dict()
    return out, rw_rep, sh_rep


def emulate(seed: int) -> dict:
    """Score every candidate of one seed exactly as vec.evaluate.score_board does, with the anchors computed once."""
    d = _DATA
    cfg, board = d["cfg"], d["board"]
    t0 = time.perf_counter()
    cands, rw_rep, sh_rep = candidates(seed)
    t_build = time.perf_counter() - t0
    target_w = subsample(d["truth_stage"], cfg["frac"], seed)
    ref_w = subsample(d["ref_stage"], cfg["frac"], seed + 1)
    truth, other = split_half(target_w, seed=0)
    floor = panel_metrics(TASK, ref_w, truth, ref_w, seed=seed, build=cfg["build"])
    ceil = panel_metrics(TASK, other, truth, ref_w, seed=seed, build=cfg["build"])
    anchors = {m: {"floor": floor.get(m), "ceiling": ceil.get(m)} for m in RANKED_METRICS[TASK]}
    anchors = {m: a for m, a in anchors.items() if a["floor"] is not None and a["ceiling"] is not None}
    rows = []
    for name, pred in cands.items():
        vals = panel_metrics(TASK, pred, truth, ref_w, seed=seed, build=cfg["build"])
        skills = skills_from_anchors(TASK, vals, anchors)
        score, groups = aggregate(TASK, skills)
        row = {"seed": seed, "candidate": name, "score": score, "n_pred": pred.n}
        row.update({f"group:{k}": v for k, v in groups.items()})
        row.update({f"skill:{m}": skills[m] for m in RANKED_METRICS[TASK]})
        row.update({m: vals.get(m) for m in (*RANKED_METRICS[TASK], *EXTRA_METRICS)})
        rows.append(row)
    check = None
    if seed == cfg["check_seed"]:  # the cached-anchor loop must equal vec.evaluate.score_board
        ref_score = score_board(board, cands["copy_last"], d["truth_stage"], d["ref_stage"], seed=seed,
                                frac=cfg["frac"], build=cfg["build"]).score
        mine = next(r["score"] for r in rows if r["candidate"] == "copy_last")
        check = {"seed": seed, "score_board": ref_score, "loop": mine, "abs_diff": abs(ref_score - mine)}
    return {"seed": seed, "rows": rows, "anchors": anchors, "truth_cells": truth.n, "ref_cells": ref_w.n,
            "reweight": rw_rep, "shift": sh_rep, "check": check, "seconds_build": t_build,
            "seconds_total": time.perf_counter() - t0}


def paired(df: pd.DataFrame, base: str) -> pd.DataFrame:
    piv = df.pivot(index="seed", columns="candidate", values="score")
    diff = piv.sub(piv[base], axis=0)
    return pd.DataFrame({f"d_{base}_mean": diff.mean(), f"d_{base}_sd": diff.std(ddof=1),
                         f"d_{base}_min": diff.min(), f"d_{base}_max": diff.max(),
                         f"wins_vs_{base}": (diff > 0).sum()})


def format_smoke(cfg: dict, out_dir: Path) -> dict:
    """build_submission + write_checked on real data (oracle shares, seed 0), compared with the reweight path; the
    file goes to runs/ and is deleted after the check."""
    d = _DATA
    board = d["board"]
    shares, within = d["share_sets"]["oracle"]
    ref = d["ref"].copy()
    sub, _ = build_submission(ref, board, shares, cfg["n_cells"], seed=0, ref_source=cfg["ref_source"],
                              within=within)
    pop, _ = reweight(d["ref"], shares, cfg["n_cells"], 0, within=within, fallback=d["fallback"])
    same = bool(np.array_equal(sub.X, pop.X.toarray()))
    path = REPO_ROOT / "runs" / "t1_scope" / f"format_smoke_{out_dir.name}.h5ad"
    check = write_checked(sub, path, board)
    size = path.stat().st_size
    path.unlink()
    stamp_path(path).unlink(missing_ok=True)
    return {"ok": check.ok, "warnings": check.warnings, "n_cells": check.n_cells, "n_genes": int(sub.n_vars),
            "dtype": str(sub.X.dtype), "file_mb": round(size / 2**20, 1), "same_rows_as_reweight": same,
            "board": board.key, "limits": [board.min_cells, board.max_cells]}


@app.command()
def main(
    ref: Path = typer.Option(Path("data/raw/T1/E8.5_RNA.h5ad"), help="reference stage (copy_last, cell source)"),
    truth_stage: Path = typer.Option(Path("data/raw/T1/E9.5_RNA.h5ad"), help="stage the truth library comes from"),
    truth_library: str = typer.Option("2", help="obs-name suffix of the truth library; the others are the sibling"),
    momentum_start: Path = typer.Option(None, help="stage before the reference, for real momentum arms"),
    contamination: float = typer.Option(0.05, help="share of the truth drawn from the contamination labels"),
    contamination_labels: str = typer.Option("Foregut,Epithelium", help="foregut-type labels (surrounding), checked "
                                             "at load time: the truth library needs cells of them, and the momentum "
                                             "endpoints (sibling libraries, reference, momentum start) --min-cells "
                                             "in all"),
    n_cells: int = typer.Option(5000, help="cells per candidate (within the T1:val limits)"),
    seeds: str = typer.Option("0,1,2,3,4,5,6,7,8,9", help="emulator seeds (>= 8)"),
    cm_delta: float = typer.Option(0.10, help="mis-specification of the CM share of heart cells"),
    blood_delta: float = typer.Option(0.10, help="blood over-prediction (share of all cells)"),
    dampings: str = typer.Option("0.5,1.0", help="shift dampings"),
    flip_dampings: str = typer.Option("0.5", help="dampings also run with the sign reversed"),
    eps: float = typer.Option(0.5, help="CP10k pseudocount of the shift's log ratios"),
    min_cells: int = typer.Option(20, help="cells per type at each stage for a shift (also the contamination labels' "
                                  "minimum at each momentum endpoint)"),
    standin_momentum: bool = typer.Option(False, help="also shift types drawn from stand-in labels at a stage (off by "
                                          "default, as in vec.t1_scope.momentum_shift)"),
    frac: float = typer.Option(SUBSAMPLE_FRAC, help="emulator working-copy fraction"),
    build: str = typer.Option(DEFAULT_BUILD, help="scorer build"),
    truth_seed: int = typer.Option(0, help="seed for the contamination draw (fixed across emulator seeds)"),
    jobs: int = typer.Option(5, help="worker processes (one seed each)"),
    board_key: str = typer.Option("T1:val", help="board for the gene order and cell limits"),
    tag: str = typer.Option("", help="results subdirectory (default pseudo_<ref>_<truth>_lib<truth-library>)"),
    smoke: bool = typer.Option(True, help="also run build_submission + the format check on real data"),
    small_n: int = typer.Option(2000, help="cell count of the oracle_n<k> ablation"),
    only: str = typer.Option("", help="comma-separated candidates to run (copy_last and oracle always run)"),
):
    """Score copy_last, chamber re-weighting and the shift arms on a chamber pseudo-task (default E8.5 -> E9.5)."""
    seed_list = [int(s) for s in seeds.split(",") if s.strip()]
    if len(seed_list) < 8:
        console.print(f"[yellow]warning[/] {len(seed_list)} seeds; the protocol asks for at least 8")
    t0 = time.perf_counter()
    git_at_start = git_state()  # the procedure files as loaded; checked again before anything is written
    spec = load_spec()
    board = get_board(board_key)
    if not board.min_cells <= n_cells <= (board.max_cells or n_cells):
        raise typer.BadParameter(f"n_cells {n_cells} outside {board_key} limits [{board.min_cells}, {board.max_cells}]")
    contam = [s.strip() for s in contamination_labels.split(",") if s.strip()]
    chamber = spec.in_scope("chamber")
    stage_of = {p: p.name.split("_")[0].removesuffix(".h5ad") for p in (ref, truth_stage, momentum_start)
                if p is not None}
    out_dir = RESULTS / (tag or f"pseudo_{stage_of[ref]}_{stage_of[truth_stage]}_lib{truth_library}")

    a = load_stage(ref, lineage=False)
    b = load_stage(truth_stage, lineage=False)
    for x in (a, b):
        if list(x.var_names) != board.genes:
            raise typer.BadParameter("stage genes differ from the board's gene order")
    ref_source, truth_source = ref.name.removesuffix(".h5ad"), truth_stage.name.removesuffix(".h5ad")
    classify(a, source=ref_source)
    classify(b, source=truth_source)
    libs = sorted(set(library_of(b.obs_names)))
    if truth_library not in libs or len(libs) < 2:
        raise typer.BadParameter(f"truth stage libraries {libs}; need at least two, including {truth_library!r}")
    sib_library = ",".join(x for x in libs if x != truth_library)
    mstart, mstart_source = None, None
    if momentum_start is not None:
        mstart = load_stage(momentum_start, lineage=False)
        if list(mstart.var_names) != board.genes:
            raise typer.BadParameter("momentum-start genes differ from the board's gene order")
        mstart_source = momentum_start.name.removesuffix(".h5ad")
        classify(mstart, source=mstart_source)
    in_truth = library_of(b.obs_names) == truth_library
    truth_lib = {f"truth {stage_of[truth_stage]} library {truth_library}": b[in_truth]}
    # momentum endpoints: the shift arms run reference -> sibling, the momentum arms momentum start -> reference
    endpoints = {f"sibling {stage_of[truth_stage]} libraries {sib_library}": b[~in_truth],
                 f"reference {stage_of[ref]}": a}
    if mstart is not None:
        endpoints[f"momentum start {stage_of[momentum_start]}"] = mstart
    # the contamination labels are unused without contamination
    contam_cells = check_contamination(truth_lib, endpoints, contam, min_cells) if contamination > 0 else None
    truth, truth_info = make_truth(b, truth_library, chamber, contamination, contam, truth_seed)
    sibling = b[~in_truth].copy()

    truth_shares = load_shares(composition(truth), spec)
    ref_shares = load_shares(composition(a), spec)
    within = {"surrounding": contam}
    cm_heart = share_of(truth_shares, "cm", spec.heart)
    blood_new = truth_shares["blood"] + blood_delta
    no_surr = {c: s for c, s in ref_shares.items() if c != "surrounding"}
    no_surr = {c: s / sum(no_surr.values()) for c, s in no_surr.items()}
    share_sets = {
        "oracle": (truth_shares, within),
        "oracle_any_surr": (truth_shares, None),
        "cm_minus": (with_share(truth_shares, "cm", max(0.0, cm_heart - cm_delta)), within),
        "cm_plus": (with_share(truth_shares, "cm", min(1.0, cm_heart + cm_delta)), within),
        "drop_surr_only": (no_surr, None),
        "blood_plus": (with_share(truth_shares, "blood", blood_new, within=list(spec.categories)), within),
        "oracle_no_epi": ({c: s / (1.0 - truth_shares["epicardium"]) if c != "epicardium" else 0.0
                           for c, s in truth_shares.items()}, within),
    }
    held = pd.read_csv(REPO_ROOT / "experiments" / "t1_direction" / "results" / "held_flat_genes.csv")["gene"]
    cfg = {"n_cells": n_cells, "frac": frac, "build": build, "dampings": [float(x) for x in dampings.split(",") if x],
           "flip_dampings": [float(x) for x in flip_dampings.split(",") if x], "eps": eps, "min_cells": min_cells,
           "check_seed": seed_list[0], "ref_source": ref_source, "small_n": small_n,
           "standin_momentum": standin_momentum, "only": [x.strip() for x in only.split(",") if x.strip()]}
    _DATA.update(ref=a, sibling=sibling, cfg=cfg, board=board, share_sets=share_sets, within=within,
                 fallback=spec.fallbacks(ref_source), sib_fallback=spec.fallbacks(truth_source),
                 momentum_start=mstart, mstart_fallback=spec.fallbacks(mstart_source),
                 held_flat=held.astype(str).tolist(), truth_stage=to_stage(truth), ref_stage=to_stage(a))
    console.print(f"loaded in {time.perf_counter() - t0:.0f} s: ref {a.n_obs} cells, truth {truth.n_obs} "
                  f"(library {truth_library}: {truth_info['n_chamber']} chamber + {truth_info['n_contamination']} "
                  f"contamination), sibling library {sib_library} {sibling.n_obs} cells; CM share of heart cells in "
                  f"the truth {cm_heart:.3f}")

    results = []
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=max(1, jobs), mp_context=ctx) as pool:
        for res in pool.map(emulate, seed_list):
            results.append(res)
            best = max(res["rows"], key=lambda r: r["score"])
            console.print(f"seed {res['seed']}: {res['seconds_total']:.0f} s; copy_last "
                          f"{next(r['score'] for r in res['rows'] if r['candidate'] == 'copy_last'):.2f}, best "
                          f"{best['candidate']} {best['score']:.2f}")

    df = pd.DataFrame([r for res in results for r in res["rows"]])
    metric_cols = [c for c in df.columns if c.startswith(("skill:", "group:")) or c in RANKED_METRICS[TASK]]
    summ = df.groupby("candidate").agg(score_mean=("score", "mean"), score_sd=("score", "std"),
                                       n_seeds=("seed", "nunique"))
    summ = summ.join(df.groupby("candidate")[metric_cols].mean())
    summ = summ.join(paired(df, "copy_last")).join(paired(df, "oracle"))
    order = list(df["candidate"].drop_duplicates())
    summ = summ.loc[order]

    comp_rows = [{"set": "truth", **truth_shares}, {"set": f"reference_{stage_of[ref]}", **ref_shares}]
    for name, (shares, _) in share_sets.items():
        comp_rows.append({"set": name, **shares})
    first = results[0]
    reuse = {name: {"reused": rep["reused"], "rows": rep["rows"]} for name, rep in first["reweight"].items()}
    anchors = pd.DataFrame([{"seed": r["seed"], "truth_cells": r["truth_cells"], "ref_cells": r["ref_cells"],
                             **{f"{m}:{k}": v for m, a in r["anchors"].items() for k, v in a.items()}}
                            for r in results])
    files = {"scores.csv": df.to_csv(index=False), "summary.csv": summ.to_csv(),
             "compositions.csv": pd.DataFrame(comp_rows).to_csv(index=False),
             "reweight_seed0.json": json.dumps(reuse, indent=2, default=float) + "\n",
             "shift_seed0.json": json.dumps(first["shift"], indent=2, default=float) + "\n",
             "anchors_computed.csv": anchors.to_csv(index=False)}

    smoke_rec = format_smoke(cfg, out_dir) if smoke else None
    run = {
        "date": time.strftime("%Y-%m-%d %H:%M:%S"), "git": git_at_start, "board": board_key,
        "inputs": {str(p): sha256_file(p) for p in (ref, truth_stage, momentum_start) if p is not None},
        "truth": {**truth_info, "categories": sorted(set(truth.obs["category"])), "chamber_scope": chamber,
                  "contamination": contamination, "contamination_labels": contam,
                  "contamination_label_cells": contam_cells, "truth_seed": truth_seed,
                  "n_cells": int(truth.n_obs), "cm_share_of_heart": cm_heart},
        "sibling_library": sib_library, "n_cells": n_cells, "seeds": seed_list, "frac": frac, "build": build,
        "cm_delta": cm_delta, "blood_delta": blood_delta, "shift": {"dampings": cfg["dampings"],
                                                                     "flip_dampings": cfg["flip_dampings"],
                                                                     "eps": eps, "min_cells": min_cells,
                                                                     "standin_momentum": standin_momentum,
                                                                     "n_held_flat": len(held)},
        "fallback_ref": {k: list(v) for k, v in spec.fallbacks(ref_source).items()},
        "momentum_start": None if momentum_start is None else str(momentum_start),
        "score_board_check": first["check"] or next((r["check"] for r in results if r["check"]), None),
        "format_smoke": smoke_rec, "seconds": round(time.perf_counter() - t0, 1),
        "seconds_per_seed": [round(r["seconds_total"], 1) for r in results],
    }
    files["run.json"] = json.dumps(run, indent=2, default=float) + "\n"
    write_results(out_dir, files, git_at_start)

    t = Table("candidate", "score", "sd", "d copy_last", "wins", "de_score", "de_dir", "mmd_u", "variogram")
    for name, r in summ.iterrows():
        t.add_row(name, f"{r['score_mean']:.2f}", f"{r['score_sd']:.2f}", f"{r['d_copy_last_mean']:+.2f}",
                  f"{int(r['wins_vs_copy_last'])}/{int(r['n_seeds'])}", f"{r['de_score']:.3f}",
                  f"{r['de_direction']:.3f}", f"{r['mmd_u']:.4f}", f"{r['variogram']:.6f}")
    console.print(t)
    console.print(f"score_board check: {run['score_board_check']}")
    console.print(f"format smoke: {smoke_rec}")
    console.print(f"wrote {out_dir} in {run['seconds']:.0f} s")


if __name__ == "__main__":
    app()
