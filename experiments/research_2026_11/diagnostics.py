"""Target-free diagnostics behind the screened research of 2026-11-05 (experiments/research_2026_11/REPORT.md).

Not tries. Neither diagnostic reads a target or computes a score. They read only:
- the reference stage the test scorer uses (T1: E10.5; T3: WT E8.75);
- our own upload files, archived in runs/archive/builds_88b057f;
- builds made here with the committed builders at the upload seed, written to a scratch folder and never uploaded.
Cell counts are never printed.

  t1-tail  Lever L1 of task a (the T1 tail decoder). de_score fills the predicted bottom n_dn with the genes of most
           negative dp = pb(file) - pb(reference). A truth down-gene needs pb(reference) >= 0.25, because log1p values
           are >= 0. So a predicted decrease of a gene below that level is a guaranteed miss whenever it reaches the
           bottom tail.
           For s1's bottom K genes, the diagnostic counts how many are infeasible, and splits each by cause:
           - scope: the chamber re-weighting alone (the same draw without the external direction);
           - direction: the rest.
  tails    t1-tail's infeasible-tail count, run on every archived upload file against its board's reference stage:
           - T1 s1 and s2 against E10.5 (32,285 genes);
           - the heart s1 and s2 against E10.5 MERFISH;
           - the embryo s1 and s2 against E7.5;
           - the four T3 files against WT E8.75.
           Each MERFISH board uses its 500-gene panel. (The completeness critic of the counter-review, wf_3eef9333-d9d.)
  t3-gate  Lever L2 of task a and lever 1 of task h. severity_slope floors when r2 <= 0.01. Over the truth's DE
           genes, r2 > 0.01 is equivalent to cos^2(dp, dt) > 0.01 + 0.99 cos^2(dp, 1). So a common, one-signed
           part of dp raises the agreement the truth must show.
           For s1, and for s1 rebuilt without its depth step and without its factor (the same seed, stream and
           count), the diagnostic reports:
           - cos^2(dp, 1) over candidate gene sets, and the gate bar it implies;
           - dp's rank correlation with pb(WT) (the expression-level pattern that depth produces);
           - the depth step's share of s1's dp.

Every JSON records the commit, whether the tracked tree was clean, and the sha256 of every input file (the N1/N2
precedent in experiments/harness/TRIES.md). Rows "diagnostic, not a try" in TRIES.md register these runs.

Run on the DGX from the repo root:
    uv run --frozen --no-sync python -I -B experiments/research_2026_11/diagnostics.py t1-tail --out DIR
    uv run --frozen --no-sync python -I -B experiments/research_2026_11/diagnostics.py t3-gate --out DIR
    uv run --frozen --no-sync python -I -B experiments/research_2026_11/diagnostics.py tails --out DIR
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import sparse, stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "t3_decomp"))

from vec import harness, t1_scope  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.prep import load_stage  # noqa: E402

ARCHIVE = REPO / "runs" / "archive" / "builds_88b057f"
SEED = harness.BUILD_SEED_OFFSET + 1
KS = (25, 50, 100, 200, 400, 800, 1600)


def pseudobulk(X) -> np.ndarray:
    """Per-gene mean of log1p CP10k over cells (the scorer's pb)."""
    X = X.tocsr() if sparse.issparse(X) else np.asarray(X)
    return np.asarray(X.mean(axis=0)).ravel().astype(np.float64)


def gene_sd(X) -> np.ndarray:
    X = X.tocsr() if sparse.issparse(X) else np.asarray(X)
    m = np.asarray(X.mean(axis=0)).ravel()
    m2 = np.asarray(X.multiply(X).mean(axis=0)).ravel() if sparse.issparse(X) else (X * X).mean(axis=0)
    return np.sqrt(np.maximum(m2 - m * m, 0.0))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(inputs: list[Path]) -> dict:
    git = ["git", "-C", str(REPO)]
    head = subprocess.run([*git, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run([*git, "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                           text=True).stdout.strip()
    return {"commit": head, "tracked_tree_clean": not dirty,
            "inputs_sha256": {str(p.relative_to(REPO)): sha256(p) for p in inputs}}


def tail_rows(dp: np.ndarray, pb_ref: np.ndarray, ks) -> list[dict]:
    order = np.argsort(dp, kind="stable")
    infeas = pb_ref < 0.25
    return [{"K": k, "infeasible": int(infeas[order[:k]].sum()), "kth_dp": float(dp[order[k - 1]])} for k in ks]


def cos2_ones(v: np.ndarray) -> float:
    n = v.size
    return float(v.sum() ** 2 / (n * float(v @ v))) if float(v @ v) > 0 else float("nan")


# ---------------------------------------------------------------------------------------------------------------------
def t1_tail(out: Path) -> dict:
    import anndata as ad

    board = load_boards()["T1:test"]
    ref = load_stage(harness.DATA_DIR / "T1/E10.5_RNA.h5ad", lineage=False)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E10.5_RNA")
    n = harness.count_for(board, ref.n_obs)  # never printed
    scope, _ = t1_scope.build_submission(ref, board, shares, n, SEED, held_flat=t1_scope.held_flat_genes(),
                                         ref_source="E10.5_RNA")
    s1 = ad.read_h5ad(ARCHIVE / "builds" / "T1_s1.h5ad")
    genes = list(s1.var_names)
    if list(scope.var_names) != genes:
        raise SystemExit("scope-only build and s1 differ in gene order")
    pos = {g: i for i, g in enumerate(ref.var_names)}
    Xr = ref.X.tocsr()[:, [pos[g] for g in genes]]
    pb_ref = pseudobulk(Xr)
    se10 = gene_sd(Xr) / np.sqrt(0.1 * ref.n_obs)  # the scorer's 10% reference copy
    dp_s1 = pseudobulk(s1.X) - pb_ref
    dp_scope = pseudobulk(scope.X) - pb_ref
    dp_dir = dp_s1 - dp_scope
    infeas = pb_ref < 0.25
    infeas_m = pb_ref < 0.25 + 2 * se10
    order = np.argsort(dp_s1, kind="stable")  # most negative first
    rows = []
    for k in KS:
        b = order[:k]
        scope_driven = np.abs(dp_scope[b]) >= 0.5 * np.abs(dp_s1[b])
        rows.append({"K": k,
                     "infeasible": int(infeas[b].sum()), "infeasible_with_margin": int(infeas_m[b].sum()),
                     "infeasible_scope_driven": int((infeas[b] & scope_driven).sum()),
                     "infeasible_direction_driven": int((infeas[b] & ~scope_driven).sum()),
                     "kth_dp": float(dp_s1[b[-1]])})
    feas_depth = {f"feasible_dp_le_{x}": int(((dp_s1 <= -x) & ~infeas).sum()) for x in (0.05, 0.1, 0.25)}
    infeas_depth = {f"infeasible_dp_le_{x}": int(((dp_s1 <= -x) & infeas).sum()) for x in (0.02, 0.05, 0.1)}
    res = {"diagnostic": "t1-tail", "reference": "T1 E10.5 (all cells, as the scorer's reference)",
           "provenance": provenance([harness.DATA_DIR / "T1/E10.5_RNA.h5ad", ARCHIVE / "builds" / "T1_s1.h5ad"]),
           "file": "runs/archive/builds_88b057f/builds/T1_s1.h5ad", "scope_only": "same draw, no direction",
           "n_genes": len(genes), "genes_feasible_down": int((~infeas).sum()), "tail": rows,
           **feas_depth, **infeas_depth,
           "spearman_dp_scope_vs_dp_dir_on_bottom_400": float(stats.spearmanr(dp_scope[order[:400]],
                                                                                dp_dir[order[:400]])[0])}
    return res


# ---------------------------------------------------------------------------------------------------------------------
def t3_gate(out: Path, threads: int) -> dict:
    import build_submission as bs
    import t3_decompose as t3d

    decision = REPO / "experiments/t3_decomp/results/gata4/gata4_decisions.json"
    run_dir = REPO / "runs/t3_decomp/gata4_20261028T122425Z"
    wt_path = harness.DATA_DIR / "T2/E8.75.h5ad"
    board = bs.board_spec("T3:bcatenin")
    n = bs.resolve_n_cells("max", board)  # never printed
    pf = bs.preflight(decision, "allowed", "s1", board, n)
    spec, dec, rules, preset = pf["spec"], pf["dec"], pf["rules"], pf["preset"]
    stage_paths = None
    if spec.get("stage"):
        stage_paths = {s: harness.DATA_DIR / "T2" / f"{s}.h5ad" for s in rules["presets"][preset]["stage_pools"]}
    ref = bs.load_reference(wt_path, preset, rules, stage_paths, threads)
    gpos = {g: i for i, g in enumerate(ref.genes)}
    panel = [g for g in board.genes]
    cols = [gpos[g] for g in panel]
    wt = load_stage(wt_path)  # as the scorer loads the reference: as released
    wpos = {g: i for i, g in enumerate(wt.var_names)}
    pb_wt = pseudobulk(wt.X.tocsr()[:, [wpos[g] for g in panel]])

    def variant(name: str, depth: bool, factor: bool) -> np.ndarray:
        sp = json.loads(json.dumps(spec))
        sp["depth"] = bool(depth) and bool(spec.get("depth"))
        kw, _ = bs.inputs_for(sp, dec, ref, run_dir)
        if not factor:
            kw["phi"], kw["w"] = None, 1.0
        built = bs.construct(ref, n, t3d.stable_rng(SEED, "submission:allowed:s1"), **kw)
        return pseudobulk(built.X[:, cols]) - pb_wt

    dps = {"s1_as_built": variant("s1", True, True), "s1_no_depth": variant("nd", False, True),
           "s1_no_factor": variant("nf", True, False), "composition_only": variant("c", False, False)}
    import anndata as ad
    arch = {}
    for name in ("allowed_s1", "refused_s1", "barred_s1", "barred_s2"):
        a = ad.read_h5ad(ARCHIVE / "t3_build" / f"{name}.h5ad")
        p = {g: i for i, g in enumerate(a.var_names)}
        arch[name] = pseudobulk(a.X[:, [p[g] for g in panel]] if sparse.issparse(a.X)
                                else np.asarray(a.X)[:, [p[g] for g in panel]]) - pb_wt
    check = float(np.corrcoef(dps["s1_as_built"], arch["allowed_s1"])[0, 1])
    feas = pb_wt >= 0.25
    rep = {}
    for name, dp in {**dps, **{f"archived_{k}": v for k, v in arch.items()}}.items():
        top = np.argsort(-np.abs(dp))[: max(1, dp.size // 5)]
        sets = {"all_genes": np.arange(dp.size), "feasible_down_or_up": np.flatnonzero(feas | (dp > 0)),
                "top20pct_abs_dp": top}
        r = {"mean_dp": float(dp.mean()), "sd_dp": float(dp.std()),
             "spearman_dp_vs_pb_wt": float(stats.spearmanr(dp, pb_wt)[0])}
        for sname, idx in sets.items():
            c2 = cos2_ones(dp[idx])
            r[f"cos2_ones_{sname}"] = c2
            r[f"gate_bar_cos2_dp_dt_{sname}"] = 0.01 + 0.99 * c2
        rep[name] = r
    d_depth = dps["s1_as_built"] - dps["s1_no_depth"]
    s1v = dps["s1_as_built"]
    res = {"diagnostic": "t3-gate", "reference": "WT E8.75 as released (the scorer's reference)",
           "provenance": provenance([wt_path, decision, *[ARCHIVE / "t3_build" / f"{n}.h5ad" for n in arch]]),
           "spec_allowed_s1": {k: spec.get(k) for k in ("depth", "stage", "coverage")},
           "rebuild_matches_archived_s1_corr": check, "variants": rep,
           "depth_step": {"share_of_s1_dp_energy": float(d_depth @ d_depth / (s1v @ s1v)),
                          "cos2_ones": cos2_ones(d_depth),
                          "spearman_vs_pb_wt": float(stats.spearmanr(d_depth, pb_wt)[0]),
                          "corr_with_s1_dp": float(np.corrcoef(d_depth, s1v)[0, 1])}}
    return res


def tails(out: Path) -> dict:
    import anndata as ad

    boards = load_boards()
    jobs = [  # (file, board key, reference file, Ks)
        ("builds/T1_s1.h5ad", "T1:test", "T1/E10.5_RNA.h5ad", KS),
        ("builds/T1_s2.h5ad", "T1:test", "T1/E10.5_RNA.h5ad", KS),
        ("builds/T2heart_s1.h5ad", "T2:heart:test_extrap", "T2/E10.5.h5ad", (20, 40, 80, 160)),
        ("builds/T2heart_s2.h5ad", "T2:heart:test_extrap", "T2/E10.5.h5ad", (20, 40, 80, 160)),
        ("builds/T2embryo_s1.h5ad", "T2:embryo:test_interp", "T2/E7.5.h5ad", (20, 40, 80, 160)),
        ("builds/T2embryo_s2.h5ad", "T2:embryo:test_interp", "T2/E7.5.h5ad", (20, 40, 80, 160)),
        *[(f"t3_build/{n}.h5ad", "T3:bcatenin", "T2/E8.75.h5ad", (20, 40, 80, 120))
          for n in ("allowed_s1", "refused_s1", "barred_s1", "barred_s2")],
    ]
    refs: dict[str, tuple[list[str], np.ndarray]] = {}
    rows, inputs = {}, []
    for rel, key, ref_rel, ks in jobs:
        board = boards[key]
        panel = list(board.genes)
        if ref_rel not in refs:
            r = load_stage(harness.DATA_DIR / ref_rel, lineage=False) if ref_rel.startswith("T1/") else \
                load_stage(harness.DATA_DIR / ref_rel)
            pos = {g: i for i, g in enumerate(r.var_names)}
            refs[ref_rel] = (panel, pseudobulk(r.X.tocsr()[:, [pos[g] for g in panel]]))
            inputs.append(harness.DATA_DIR / ref_rel)
        rpanel, pb_ref = refs[ref_rel]
        if rpanel != panel:
            raise SystemExit(f"{rel}: panel differs from the cached reference's")
        a = ad.read_h5ad(ARCHIVE / rel)
        p = {g: i for i, g in enumerate(a.var_names)}
        X = a.X[:, [p[g] for g in panel]]
        dp = pseudobulk(X) - pb_ref
        rows[rel] = {"board": key, "reference": ref_rel, "genes": len(panel),
                     "genes_feasible_down": int((pb_ref >= 0.25).sum()), "tail": tail_rows(dp, pb_ref, ks)}
        inputs.append(ARCHIVE / rel)
    return {"diagnostic": "tails", "provenance": provenance(inputs), "files": rows}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("which", choices=("t1-tail", "t3-gate", "tails"))
    ap.add_argument("--out", required=True, help="folder for the JSON result")
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {"t1-tail": lambda: t1_tail(out), "t3-gate": lambda: t3_gate(out, a.threads),
           "tails": lambda: tails(out)}[a.which]()
    path = out / f"{a.which}.json"
    path.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
