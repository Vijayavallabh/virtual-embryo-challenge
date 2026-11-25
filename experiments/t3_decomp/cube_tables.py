#!/usr/bin/env python
"""Row R72's sign cube (the captain's order of 2026-11-25; T3_MOON5_design.md with the fixes F1 to F10 of
T3_MOON5_challenge.md and T3_MOON5_check2.md): its four WT-only tables and the target-free gate G0 (a to g).

  tables --data DATA --out DIR
      From released files only (WT E8.75, E9.5 and E8.0; the released Mab21l2 knockout and the two unlocked Gata4
      replicates for the two pinned tables' label masks; rules.yaml), in the frozen run's context as
      t3_decompose.build_context makes it (the gata4 preset's WT E8.75 pool, label mode, vocabulary and harmoniser):
      WT E9.5 and WT E8.0 are called into the vocabulary as t3_decompose.transfer_factors calls a knockout file, and
        S_M, Y_M = gene_factor(WT E9.5, WT E8.0 over WT E8.75) on transfer_mab21l2_all's mask (labels with >= 100
                   Mab21l2 cells and >= 100 WT cells), so on that table's own den
        S_G, Y_G = the same on gata4_pooled's mask (strict non-Mesp1 labels with >= 20 cells in the pooled replicates
                   and in WT), so on that table's own den
      with the same weights (WT shares), eps and clip as the pinned tables; a masked label whose stage sample holds
      under rules.yaml's stage min_cells (30) takes WT E8.75's own label mean as its numerator, so the den is
      unchanged. Writes DIR/<table>.csv (gene, factor; WT E8.75's gene order, as the frozen factor files) and
      DIR/MANIFEST.json (each table's sha256, the masks, the inputs' sha256, the code's sha256, commit, clean, machine).
      No pathway, gene set, model or outside sample enters; no target is compared with anything.
  g0 --data DATA --run-dir RUN_DIR --decision DEC --manifest MANIFEST --manifest-sha256 SHA --out DIR
      G0, target-free (WT files, the committed tables, C2's drawn cells in memory; no truth, no score; no file but its
      reading): DIR/g0.json and DIR/g0.txt, whose first line is the one fixed outcome line. Exit 0 CONTINUE, 5 STOP,
      2 refused (a dirty tree, an output already there, the manifest or a table off its sha256), 1 a crash (masked).
      G0a  transfer_mab21l2_all and gata4_pooled recomputed from released files on the masks this context identifies:
           the masks equal the manifest's and max |ln recomputed - ln pinned| <= 1e-6 for each, else STOP (the dens
           are not pinned); prints the E8.0 and E9.5 cells called into each masked label
      G0b  A + X_S + X_M + X_G = c (l_M + l_G) to 1e-12, and against ln phi_C2; the genes at a clip in any table; at
           ppp the vertex path gives |k - 1| <= 1e-9 (then k = 1) and C2's factor array exactly (the file-level
           identity is the launcher's rebuild of C2's bytes 47111ac1...); else STOP
      G0c  SD(Q(X_j)) / SD(ln phi_C2) >= 0.15 for each axis, else that axis is held at '+'
      G0d  every vertex's k_v, range, den total and each held label's expected library against C2's: the builder's
           guards (k_v in [0.5, 2], the hold converged) and every label with a cell in C2's draw (row R72's every
           label) within +-1% of C2's library, else the vertex is not built
      G0e  per axis, the rank correlation (Spearman) over the board's panel of dp(C2 with the axis alone reflected) -
           dp(C2) against the dp change of a uniform depth change of C2's cells to WT E9.5's and to WT E8.0's median
           depth over WT E8.75's: |rho| <= 0.3 on both, else held at '+' (depth by another route); Pearson reported
      G0f  the same against the dp change of row R36's composition block on C2's command and of C2's cells redrawn
           at WT E8.0's and WT E9.5's label shares: |rho| <= 0.3 on all three, else held at '+' (composition)
      G0g  on l_M's den, s(x) = <ln x - ln WT8.0, ln WT9.5 - ln WT8.0> / |ln WT9.5 - ln WT8.0|^2: 0 < s(E8.25_late)
           < s(E8.75) < 1, |s(E8.75) - 0.5| <= 0.25, E8.5 between them; else X_S is held at '+' (N9's known lag)
      STOP if G0a or G0b fails or neither S nor M is left; otherwise CONTINUE with the axes left and the vertices
      G0 leaves (those reflecting only axes left, and passing G0d).
No log holds a score or a cell count of a test board: C2's drawn cells are counted by code only.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import re
import socket
import subprocess
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in (str(HERE), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import sign_cube as sc  # noqa: E402

SCRIPT = "experiments/t3_decomp/cube_tables.py"
SCRIPT_PATH = Path(__file__).resolve()
TABLES = ("S_M", "Y_M", "S_G", "Y_G")
STAGE_OF = {"S": "E9.5", "Y": "E8.0"}  # the older and the younger pole
MASK_TABLE = {"M": "transfer_mab21l2_all", "G": "gata4_pooled"}
WT_FILES = {"E8.75": "T2/E8.75.h5ad", "E9.5": "T2/E9.5.h5ad", "E8.0": "T2/E8.0.h5ad",
            "E8.25_late": "T2/E8.25_late.h5ad", "E8.5": "T2/E8.5.h5ad"}
KO_FILES = {"mab21l2": "T3/E9.5_mab21l2_ko.h5ad", "rep1": "T3/E8.75_gata4_ko_1.h5ad",
            "rep2": "T3/E8.75_gata4_ko_2.h5ad"}
PRESET = "gata4"
TRANSFER_MIN_CELLS = 100  # t3_decompose.transfer_factors: both sides >= 100 cells
MANIFEST = "MANIFEST.json"
DEFAULT_OUT = HERE / "results" / "t3_cube"
STREAM, SEED, BRANCH, SLOT = "submission:barred:s2", 1_000_001, "barred", "s2"  # C2's draw (the barred s2)
G0A_TOL, G0B_TOL = 1e-6, 1e-12
G0G_AXIS = ("E8.25_late", "E8.5", "E8.75")
EXIT_STOP, EXIT_REFUSED, EXIT_CRASH = 5, 2, 1
OUTCOME_STOP = "Cube G0: STOP ({why}): nothing is read or built."
OUTCOME_GO = ("Cube G0: CONTINUE: axes left {kept}; held at '+' {held}; vertices G0d refuses {refused}; vertices G0 "
              "leaves {leaves}.")


class Refused(Exception):
    pass


def _bs():
    import build_submission as bs

    return bs


def _btk():
    import build_two_knockout as btk

    return btk


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _git(*args: str) -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, timeout=60)
    except Exception:  # noqa: BLE001
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def stamp() -> dict:
    """The commit, the clean flag (tracked files: git status --porcelain --untracked-files=no prints nothing), the
    machine and the time."""
    st = _git("status", "--porcelain", "--untracked-files=no")
    return {"commit": _git("rev-parse", "HEAD"), "tree_clean": st == "", "machine": socket.gethostname(),
            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}


_DECIMAL = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[eE][-+]?\d+)?|[-+]?\d+[eE][-+]?\d+|\b(?:nan|inf)\b")


def masked(text_: str) -> str:
    """A crash's text with every decimal number masked (k18_built_files.masked's rule)."""
    return _DECIMAL.sub("<num>", text_)


# ---------------------------------------------------------------------------------------------------------------------
# The frozen run's context (t3_decompose.build_context's first part, no target, no stage pool)


def wt_context(data: Path, rules: dict, n_jobs: int) -> SimpleNamespace:
    """WT E8.75 as build_context prepares it for the gata4 preset: load_stage with duplicates dropped, the preset's
    label mode, the harmoniser (rules.yaml's settings), the vocabulary from its depth-free model, the per-label mean
    within-cell proportions (wt.mu) and cell counts (wt.count)."""
    from vec.prep import load_stage

    t3d = _bs().t3d
    preset, hp = rules["presets"][PRESET], rules["harmonisation"]
    wt_ad = load_stage(Path(data) / WT_FILES["E8.75"], drop_duplicates=True)
    t3d.add_label_column(wt_ad, preset["label_mode"])
    genes = list(wt_ad.var_names)
    harm = t3d.Harmoniser(wt_ad, "t3_label", t3d.label_groups(wt_ad, "t3_label"), hp, n_jobs)
    V = t3d.make_vocab(wt_ad, "t3_label", harm.model(None), rules)
    pos = {t: i for i, t in enumerate(V.labels)}
    wlab = np.array([pos[t] for t in wt_ad.obs["t3_label"].astype(str)])
    Cw = t3d.dense_counts(wt_ad)
    return SimpleNamespace(rules=rules, hp=hp, preset=preset, genes=genes, harm=harm, V=V,
                           mu=t3d.props_type_means(Cw, wlab, V.L), count=np.bincount(wlab, minlength=V.L))


def called(ctx, path: Path, littermates: bool = False) -> tuple[np.ndarray, np.ndarray, object]:
    """A released file called into the vocabulary as transfer_factors calls a knockout file (load_stage on the WT
    genes, duplicates dropped; the harmoniser at the file's depth; no confidence filter). littermates: as
    build_context prepares a target (littermate wild-type cells split off, the preset's label column). Returns (the
    label index per cell, -1 for none; dense counts; the harmoniser's depth key)."""
    from vec.prep import load_stage

    t3d = _bs().t3d
    a = load_stage(Path(path), genes=ctx.genes, drop_duplicates=True)
    if littermates:
        lm, _, _ = t3d.littermate_split(a)
        if lm.any():
            a = a[~lm].copy()
        t3d.add_label_column(a, ctx.preset["label_mode"])
    calls, key = ctx.harm.calls(a, seed=ctx.hp["seed"])
    return t3d.calls_to_index(ctx.V, calls), t3d.dense_counts(a), key


def means(ctx, lab, C) -> tuple[np.ndarray, np.ndarray]:
    t3d = _bs().t3d
    return t3d.props_type_means(C, lab, ctx.V.L), np.bincount(lab[lab >= 0], minlength=ctx.V.L)


def ko_masks(ctx, data: Path) -> tuple[dict, dict, dict]:
    """The two pinned tables' label masks, as the frozen run made them: transfer_mab21l2_all, labels with >= 100
    Mab21l2 cells and >= 100 WT cells (transfer_factors); gata4_pooled, strict non-Mesp1 labels with >= min_cells (20)
    cells in the two replicates pooled and in WT (submission_parameters). Returns (masks, the knockouts' label means,
    the depth keys)."""
    fc = ctx.rules["components"]["factor"]
    lab, C, key_m = called(ctx, Path(data) / KO_FILES["mab21l2"])
    mu_m, cnt_m = means(ctx, lab, C)
    both = (cnt_m >= TRANSFER_MIN_CELLS) & (ctx.count >= TRANSFER_MIN_CELLS)
    labs, Cs, keys = [], [], {"mab21l2": key_m}
    for r in ("rep1", "rep2"):
        lab_r, C_r, keys[r] = called(ctx, Path(data) / KO_FILES[r], littermates=True)
        labs.append(lab_r)
        Cs.append(C_r)
    mu_g, cnt_g = means(ctx, np.concatenate(labs), np.vstack(Cs))
    ok = (cnt_g >= fc["min_cells"]) & (ctx.count >= fc["min_cells"])
    return {"M": both, "G": ok & np.asarray(ctx.V.nm, dtype=bool)}, {"M": mu_m, "G": mu_g}, keys


def stage_table(ctx, mu_x, cnt_x, mask) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """gene_factor(stage x over WT E8.75) on `mask` with the pinned tables' weights, eps and clip; a masked label with
    fewer than rules.yaml's stage min_cells cells of x takes WT E8.75's own label mean as its numerator. Returns
    (phi, den, the labels that took the reference's mean)."""
    t3d = _bs().t3d
    fc = ctx.rules["components"]["factor"]
    mn = ctx.rules["components"]["stage"]["min_cells"]
    mask = np.asarray(mask, dtype=bool)
    low = mask & (np.asarray(cnt_x) < mn)
    num = np.array(mu_x, dtype=np.float64, copy=True)
    num[low] = ctx.mu[low]
    phi, den = t3d.gene_factor(num, ctx.mu, ctx.V.p, mask, fc["eps"], tuple(fc["clip"]))
    return phi, den, low


def stage_means(ctx, data: Path, stages) -> tuple[dict, dict, dict]:
    out_mu, out_cnt, keys = {}, {}, {}
    for s in stages:
        lab, C, keys[s] = called(ctx, Path(data) / WT_FILES[s])
        out_mu[s], out_cnt[s] = means(ctx, lab, C)
    return out_mu, out_cnt, keys


def labels_of(ctx, mask) -> list[str]:
    return [t for t, b in zip(ctx.V.labels, np.asarray(mask, dtype=bool), strict=True) if b]


def write_table(path: Path, genes: list[str], phi) -> str:
    import pandas as pd

    pd.DataFrame({"gene": genes, "factor": np.asarray(phi, dtype=np.float64)}).to_csv(path, index=False)
    return sha256(path)


def read_table(path: Path, genes: list[str]) -> np.ndarray:
    """A cube table (gene, factor) in `genes`' order (refused unless the file's genes are exactly `genes`)."""
    import pandas as pd

    df = pd.read_csv(path, float_precision="round_trip")
    if list(df["gene"].astype(str)) != list(genes):
        raise Refused(f"{Path(path).name}: genes differ from WT E8.75's (order included)")
    phi = df["factor"].to_numpy(np.float64)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused(f"{Path(path).name}: values must be finite and positive")
    return phi


def load_manifest(path: Path, want_sha256: str | None = None) -> tuple[dict, str]:
    """MANIFEST.json and its sha256; refused when the sha256 is not `want_sha256` (when given) or a table's bytes are
    not the manifest's."""
    path = Path(path)
    if not path.is_file():
        raise Refused(f"{path}: no manifest")
    got = sha256(path)
    if want_sha256 is not None and got != want_sha256:
        raise Refused(f"the manifest's sha256 {got[:16]} is not {want_sha256[:16]}")
    man = json.loads(path.read_text())
    if sorted(man.get("tables", {})) != sorted(TABLES):
        raise Refused("the manifest does not hold the four tables")
    for t, rec in man["tables"].items():
        f = path.parent / rec["file"]
        if not f.is_file() or sha256(f) != rec["sha256"]:
            raise Refused(f"table {t} is not the manifest's bytes")
    return man, got


def load_tables(manifest: Path, genes: list[str], want_sha256: str | None = None) -> tuple[dict, dict, str]:
    """The four tables in `genes`' order, after load_manifest's checks; the manifest's genes must be `genes`."""
    man, got = load_manifest(manifest, want_sha256)
    if man.get("genes_sha256") != hashlib.sha256("\n".join(genes).encode()).hexdigest():
        raise Refused("the manifest's genes are not WT E8.75's")
    return {t: read_table(Path(manifest).parent / man["tables"][t]["file"], genes) for t in TABLES}, man, got


# ---------------------------------------------------------------------------------------------------------------------
# tables


def make_tables(data: Path, out: Path, threads: int = 8) -> dict:
    bs = _bs()
    st = stamp()
    if not st["tree_clean"] or not st["commit"]:
        raise Refused("the tracked tree is not clean (or git cannot tell): commit first")
    out = Path(out)
    if (out / MANIFEST).exists() or any((out / f"{t}.csv").exists() for t in TABLES):
        raise Refused(f"{out} already holds the tables (never overwritten)")
    rules = bs.load_rules()
    if rules["components"]["stage"]["min_cells"] != 30:
        raise Refused("rules.yaml's stage min_cells is not 30")
    with bs.thread_limits(threads):
        ctx = wt_context(Path(data), rules, threads)
        masks, _, ko_keys = ko_masks(ctx, Path(data))
        mu, cnt, keys = stage_means(ctx, Path(data), ("E9.5", "E8.0"))
    tabs, recs, dens = {}, {}, {}
    for t in TABLES:
        pole, m = t.split("_")
        phi, den, low = stage_table(ctx, mu[STAGE_OF[pole]], cnt[STAGE_OF[pole]], masks[m])
        tabs[t], dens[m] = phi, den
        recs[t] = {"stage": STAGE_OF[pole], "den_of": MASK_TABLE[m], "fallback_labels": labels_of(ctx, low),
                   "den_sha256": sc.factor_sha256(den), "genes_at_clip": int(((phi <= 0.25) | (phi >= 4.0)).sum())}
    out.mkdir(parents=True, exist_ok=True)
    for t in TABLES:
        recs[t] = {"file": f"{t}.csv", "sha256": write_table(out / f"{t}.csv", ctx.genes, tabs[t]), **recs[t]}
    files = {**{k: WT_FILES[k] for k in ("E8.75", "E9.5", "E8.0")}, **KO_FILES}
    man = {
        "what": "the T3 sign cube's WT-only tables (row R72; T3_MOON5_design.md 2.2 with fixes F1 to F10): S_M, Y_M = "
                "WT E9.5, WT E8.0 over WT E8.75 on transfer_mab21l2_all's den; S_G, Y_G = the same on gata4_pooled's "
                "den; t3_decompose.gene_factor (WT-share weights, eps, clip [0.25, 4]), the factor column is phi",
        "script": SCRIPT, "script_sha256": sha256(SCRIPT_PATH),
        "code_sha256": {"sign_cube.py": sha256(HERE / "sign_cube.py"), "t3_decompose.py": sha256(bs.T3D_PATH),
                        "rules.yaml": sha256(bs.RULES_PATH), "build_submission.py": sha256(bs.BUILDER_PATH)},
        "provenance": st,
        "inputs": {k: {"path": f"data/raw/{v}", "sha256": sha256(Path(data) / v)} for k, v in files.items()},
        "genes_sha256": hashlib.sha256("\n".join(ctx.genes).encode()).hexdigest(), "n_genes": len(ctx.genes),
        "vocabulary_sha256": hashlib.sha256("\n".join(ctx.V.labels).encode()).hexdigest(),
        "tables": recs,
        "masks": {m: {"table": MASK_TABLE[m], "labels": labels_of(ctx, masks[m]),
                      "den_sha256": sc.factor_sha256(dens[m])} for m in ("M", "G")},
        "mask_rules": {"M": "labels with >= 100 Mab21l2 cells and >= 100 WT E8.75 cells (transfer_factors)",
                       "G": "strict non-Mesp1 labels with >= 20 cells in the two Gata4 replicates pooled and in WT "
                            "E8.75 (submission_parameters; littermates split off)"},
        "stage_cells_in_masked_labels": {s: {m: {t: int(cnt[s][i]) for i, t in enumerate(ctx.V.labels) if masks[m][i]}
                                             for m in ("M", "G")} for s in ("E9.5", "E8.0")},
        "harmoniser_depth_keys": {**{k: (None if v is None else float(v)) for k, v in keys.items()},
                                  **{k: (None if v is None else float(v)) for k, v in ko_keys.items()}},
        "constants": {"eps": rules["components"]["factor"]["eps"], "clip": rules["components"]["factor"]["clip"],
                      "stage_min_cells": rules["components"]["stage"]["min_cells"],
                      "transfer_min_cells": TRANSFER_MIN_CELLS, "factor_min_cells":
                          rules["components"]["factor"]["min_cells"]},
    }
    (out / MANIFEST).write_text(json.dumps(man, indent=1) + "\n")
    return {"manifest": str(out / MANIFEST), "sha256": sha256(out / MANIFEST)}


# ---------------------------------------------------------------------------------------------------------------------
# G0


def spearman(a, b) -> float:
    from scipy.stats import rankdata

    ra, rb = rankdata(np.asarray(a, dtype=np.float64)), rankdata(np.asarray(b, dtype=np.float64))
    return float(np.corrcoef(ra, rb)[0, 1])


def pearson(a, b) -> float:
    return float(np.corrcoef(np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64))[0, 1])


@contextlib.contextmanager
def lognorm_spy(seen: list):
    """t3_decompose.lognorm, for this process only: the counts it receives are recorded (C2's final counts)."""
    t3d = _bs().t3d
    original = t3d.lognorm

    def spy(C):
        seen.append(np.array(C, copy=True))
        return original(C)

    t3d.lognorm = spy
    try:
        yield
    finally:
        t3d.lognorm = original


def shares(ctx_V, cnt) -> np.ndarray:
    cnt = np.asarray(cnt, dtype=np.float64)
    if not cnt.sum() > 0:
        raise Refused("a stage file has no cell called into the vocabulary")
    return cnt / cnt.sum()


def g0_outcome(res: dict) -> tuple[str, str, list[str], list[str]]:
    """The outcome from the checks: (STOP or CONTINUE, the fixed line, the axes left, the vertices G0 leaves)."""
    if not res["G0a"]["passes"]:
        return "STOP", OUTCOME_STOP.format(why="G0a: the dens are not pinned"), [], []
    if not res["G0b"]["passes"]:
        return "STOP", OUTCOME_STOP.format(why="G0b: the identity does not hold"), [], []
    held = {}
    for j in sc.AXES:
        why = [k for k in ("G0c", "G0e", "G0f") if not res[k]["axes"][j]["passes"]]
        if j == "S" and not res["G0g"]["passes"]:
            why.append("G0g")
        if why:
            held[j] = why
    kept = [j for j in sc.AXES if j not in held]
    if not ({"S", "M"} & set(kept)):
        return "STOP", OUTCOME_STOP.format(why="G0: neither S nor M is left"), kept, []
    refused = [v for v in sc.NAMED if not res["G0d"]["vertices"][v]["passes"]]
    leaves = [v for v in sc.kept_vertices(kept) if v not in refused]
    line = OUTCOME_GO.format(kept=", ".join(kept), held=(", ".join(f"{j} ({'+'.join(w)})" for j, w in held.items())
                                                         or "none"),
                             refused=", ".join(refused) or "none", leaves=", ".join(leaves) or "none")
    return "CONTINUE", line, kept, leaves


def run_g0(data: Path, run_dir: Path, decision: Path, manifest: Path, manifest_sha: str, board_key: str,
           threads: int = 16, boards_dir: Path | None = None) -> dict:
    """G0 (a to g); returns the reading (the outcome included). Writes nothing."""
    bs, btk = _bs(), _btk()
    t3d = bs.t3d
    from vec.harmonise import median_depth
    from vec.prep import load_stage

    btk.check_frozen()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    wt = Path(data) / WT_FILES["E8.75"]
    with bs.thread_limits(threads):
        ctx = wt_context(Path(data), rules, threads)
        genes = ctx.genes
        tabs, man, man_sha = load_tables(Path(manifest), genes, manifest_sha)
        # G0a
        masks, mus, ko_keys = ko_masks(ctx, Path(data))
        phi_m, rec_m = bs.load_factor(Path(run_dir), btk.TABLE_M, genes, dec)
        phi_g, rec_g = bs.load_factor(Path(run_dir), btk.TABLE_G, genes, dec)
        for nm, rec in ((btk.TABLE_M, rec_m), (btk.TABLE_G, rec_g)):
            if rec["sha256"] != btk.TABLE_SHA256[nm]:
                raise Refused(f"{nm} is not the pinned table")
        fc = rules["components"]["factor"]
        diffs, same = {}, {}
        for m, pinned in (("M", phi_m), ("G", phi_g)):
            re_, _ = t3d.gene_factor(mus[m], ctx.mu, ctx.V.p, masks[m], fc["eps"], tuple(fc["clip"]))
            diffs[m] = float(np.abs(np.log(re_) - np.log(pinned)).max())
            same[m] = labels_of(ctx, masks[m]) == man["masks"][m]["labels"]
        st_mu, st_cnt, st_keys = stage_means(ctx, Path(data), ("E9.5", "E8.0", "E8.25_late", "E8.5"))
        g0a = {"max_abs_ln_diff": {MASK_TABLE[m]: diffs[m] for m in diffs}, "tolerance": G0A_TOL,
               "masks_equal_manifest": {MASK_TABLE[m]: bool(same[m]) for m in same},
               "stage_cells_in_masked_labels": {s: {m: {t: int(st_cnt[s][i]) for i, t in enumerate(ctx.V.labels)
                                                        if masks[m][i]} for m in ("M", "G")} for s in ("E9.5", "E8.0")},
               "passes": bool(all(d <= G0A_TOL for d in diffs.values()) and all(same.values()))}
        # the axes
        l_m, l_g = np.log(phi_m), np.log(phi_g)
        phi_c2, crec = btk.combine("c2", 1.0, phi_m, phi_g)
        c = float(crec["c"])
        lt = {t: np.log(tabs[t]) for t in TABLES}
        ax = sc.cube_axes(l_m, l_g, lt["S_M"], lt["Y_M"], lt["S_G"], lt["Y_G"], c)
        tot = ax["A"] + ax["S"] + ax["M"] + ax["G"]
        # C2's draw (rows as drawn, before the factor step)
        board = bs.board_spec(board_key, boards_dir)
        n = bs.resolve_n_cells("max", board)
        got: dict = {}
        with btk.combined_loader("c2", 1.0, got):
            pf = bs.preflight(Path(decision), BRANCH, SLOT, board, n)
            ref = bs.load_reference(wt, pf["preset"], rules, None, threads)
            bs.check_inputs_against_run(bs.run_records(Path(run_dir), pf["preset"]), pf["preset"],
                                        ref.record["sha256"], {})
            kw, _ = bs.inputs_for(pf["spec"], pf["dec"], ref, Path(run_dir))
        if list(ref.genes) != genes or not np.array_equal(kw["phi"], phi_c2) or kw["w"] != 1.0:
            raise Refused("C2's build does not apply C2's combined factor over WT E8.75's genes")
        if list(ref.V.labels) != list(ctx.V.labels):
            raise Refused("the build's vocabulary is not the context's (label order)")
        rng = lambda name: t3d.stable_rng(SEED, name)  # noqa: E731
        drawn = bs.construct(ref, n, rng(STREAM), **{**kw, "phi": None})
        held, M = sc.label_rows(ref.pool.C[drawn.row], drawn.lab)
        B = sc.basis(M, phi_c2)
        Q = sc.projector(B)
        q = {j: Q(ax[j]) for j in sc.AXES}
        _, den = _wt_den(wt, genes)
        # G0b
        clip_counts = {t: int(((np.exp(v) <= 0.25 + 1e-12) | (np.exp(v) >= 4.0 - 1e-12)).sum())
                       for t, v in {"l_M": l_m, "l_G": l_g, **lt}.items()}
        id_err = float(np.abs(tot - c * (l_m + l_g)).max())
        id_err_c2 = float(np.abs(tot - np.log(phi_c2)).max())
        try:
            phi_ppp, rppp = sc.vertex_factor(phi_c2, q, sc.IDENTITY, M, den, B)
            ppp_k, ppp_same = rppp["k"], bool(np.array_equal(phi_ppp, phi_c2))
        except sc.CubeRefused as e:
            ppp_k, ppp_same, rppp = None, False, {"refused": str(e)}
        g0b = {"max_abs_identity_error": id_err, "max_abs_error_against_ln_phi_c2": id_err_c2, "tolerance": G0B_TOL,
               "genes_at_a_clip": clip_counts, "ppp_k": ppp_k, "ppp_is_c2_factor": ppp_same,
               "passes": bool(id_err <= G0B_TOL and ppp_k == 1.0 and ppp_same)}
        # G0c
        sd_c2 = sc.sd(np.log(phi_c2))
        g0c = {"bar": 0.15, "axes": {j: {"sd_q_over_sd_ln_phi_c2": sc.sd(q[j]) / sd_c2,
                                         "passes": bool(sc.sd(q[j]) / sd_c2 >= 0.15)} for j in sc.AXES}}
        # G0d
        verts, facs = {}, {}
        for v in sc.VERTICES:
            try:
                facs[v], r = sc.vertex_factor(phi_c2, q, v, M, den, B)
                verts[v] = {k: r[k] for k in ("k", "min", "max", "den_dev", "max_abs_library_dev", "clipped_low",
                                              "clipped_high", "sd_ln_factor", "factor_sha256")} | {
                    "hold_rounds": r["hold"]["rounds"], "passes": sc.within_band(r)}
            except sc.CubeRefused as e:
                verts[v] = {"refused": str(e), "passes": False}
        g0d = {"band": sc.LIB_BAND, "held_labels": int(len(held)), "vertices": verts}
        # G0e and G0f: dp changes over the board's panel (pseudobulks of log1p CP10k; the reference cancels)
        gpos = {g: i for i, g in enumerate(genes)}
        cols = [gpos[g] for g in board.genes]
        seen: list = []
        with lognorm_spy(seen):
            X_c2 = bs.construct(ref, n, rng(STREAM), **kw).X
        if len(seen) != 1:
            raise RuntimeError("C2's draw did not pass through lognorm once")
        pb_c2 = X_c2[:, cols].astype(np.float64).mean(0)
        C_fin, P_rows = seen[0], ref.pool.P[drawn.row]
        del X_c2

        def pb_of(X) -> np.ndarray:
            return np.asarray(X, dtype=np.float64)[:, cols].mean(0)

        depth = {}
        wt_depth = median_depth(load_stage(wt, drop_duplicates=True))
        for s in ("E9.5", "E8.0"):
            depth[s] = median_depth(load_stage(Path(data) / WT_FILES[s], drop_duplicates=True)) / wt_depth
        d_depth = {}
        for s, r in depth.items():
            Ct = t3d.rescale_depth(C_fin, np.full(C_fin.shape[0], r), P_rows, rng(f"cube:g0e:depth:{s}"))
            d_depth[s] = pb_of(t3d.lognorm(Ct)) - pb_c2
        block, _ = btk.composition_block(dec, rules, Path(decision), "allowed-s1")
        f_cd7, _ = bs.composition(ref.V, btk.with_composition(pf["spec"], block), {})
        d_comp = {"R36_composition": pb_of(bs.construct(ref, n, rng(STREAM), **{**kw, "f": f_cd7}).X) - pb_c2}
        for s in ("E8.0", "E9.5"):
            d_comp[f"{s}_label_shares"] = pb_of(bs.construct(ref, n, rng(STREAM),
                                                             **{**kw, "f": shares(ref.V, st_cnt[s])}).X) - pb_c2
        single = {"S": "npp", "M": "pnp", "G": "ppn"}
        g0e, g0f = {"bar": 0.3, "depth_ratios": depth, "axes": {}}, {"bar": 0.3, "axes": {}}
        for j, v in single.items():
            if v not in facs:
                g0e["axes"][j] = {"passes": False, "why": f"vertex {v} not made"}
                g0f["axes"][j] = {"passes": False, "why": f"vertex {v} not made"}
                continue
            d_j = pb_of(bs.construct(ref, n, rng(STREAM), **{**kw, "phi": facs[v]}).X) - pb_c2
            e = {s: {"spearman": spearman(d_j, d), "pearson": pearson(d_j, d)} for s, d in d_depth.items()}
            f = {s: {"spearman": spearman(d_j, d), "pearson": pearson(d_j, d)} for s, d in d_comp.items()}
            g0e["axes"][j] = {"directions": e, "passes": bool(all(abs(x["spearman"]) <= 0.3 for x in e.values()))}
            g0f["axes"][j] = {"directions": f, "passes": bool(all(abs(x["spearman"]) <= 0.3 for x in f.values()))}
        # G0g
        pos_tabs = {"E9.5": lt["S_M"], "E8.0": lt["Y_M"], "E8.75": np.zeros_like(l_m)}
        for s in ("E8.25_late", "E8.5"):
            phi_s, _, _ = stage_table(ctx, st_mu[s], st_cnt[s], masks["M"])
            pos_tabs[s] = np.log(phi_s)
        pos = sc.stage_positions(pos_tabs)
        g0g = {"positions": {s: pos[s] for s in G0G_AXIS}, "passes": sc.g0g_passes(pos)}
    res = {"G0a": g0a, "G0b": g0b, "G0c": g0c, "G0d": g0d, "G0e": g0e, "G0f": g0f, "G0g": g0g, "c": c,
           "manifest_sha256": man_sha, "harmoniser_depth_keys": {**{k: _key(v) for k, v in ko_keys.items()},
                                                                 **{k: _key(v) for k, v in st_keys.items()}},
           "inputs_sha256": {"decision": sha256(decision), btk.TABLE_M: rec_m["sha256"], btk.TABLE_G: rec_g["sha256"]}}
    res["outcome"], res["line"], res["axes_left"], res["vertices_left"] = g0_outcome(res)
    return res


def _key(v):
    return None if v is None else float(v)


def _wt_den(wt: Path, genes: list[str]) -> tuple[list[str], np.ndarray]:
    """build_outside_form.wt_den on WT E8.75, checked to be over `genes` in order."""
    import build_outside_form

    g, den = build_outside_form.wt_den(Path(wt))
    if list(g) != list(genes):
        raise Refused("wt_den's genes are not WT E8.75's")
    return g, np.asarray(den, dtype=np.float64)


def g0_text(res: dict) -> str:
    L = [res["line"], f"manifest sha256 {res['manifest_sha256']}; c = {res['c']:.6f}"]
    a = res["G0a"]
    L.append("G0a: max |ln recomputed - ln pinned| " + ", ".join(
        f"{k} {v:.3g}" for k, v in a["max_abs_ln_diff"].items())
             + f" (bar {a['tolerance']:g}); masks equal the manifest's " + ", ".join(
                 f"{k} {v}" for k, v in a["masks_equal_manifest"].items()) + f": {'pass' if a['passes'] else 'FAIL'}")
    for s, per in a["stage_cells_in_masked_labels"].items():
        for m, d in per.items():
            L.append(f"  {s} cells in {MASK_TABLE[m]}'s labels: " + ", ".join(f"{t} {x}" for t, x in d.items()))
    b = res["G0b"]
    L.append(f"G0b: identity error {b['max_abs_identity_error']:.3g} (bar {b['tolerance']:g}), against ln phi_C2 "
             f"{b['max_abs_error_against_ln_phi_c2']:.3g}; ppp k {b['ppp_k']}, C2's factor exactly "
             f"{b['ppp_is_c2_factor']}; "
             "genes at a clip " + ", ".join(f"{k} {v}" for k, v in b["genes_at_a_clip"].items())
             + f": {'pass' if b['passes'] else 'FAIL'}")
    L.append("G0c: SD(Q(X_j)) / SD(ln phi_C2) " + ", ".join(
        f"{j} {v['sd_q_over_sd_ln_phi_c2']:.3f} ({'pass' if v['passes'] else 'held'})"
        for j, v in res["G0c"]["axes"].items()) + " (bar 0.15)")
    L.append(f"G0d: {res['G0d']['held_labels']} labels held (every label of C2's drawn cells); band {sc.LIB_BAND:g}")
    for v, r in res["G0d"]["vertices"].items():
        if "refused" in r:
            L.append(f"  {v}: not made ({r['refused']})")
        else:
            L.append(f"  {v}: k {r['k']:.4f}, range [{r['min']:.4g}, {r['max']:.4g}], den {r['den_dev']:+.2e}, largest "
                     f"label library {r['max_abs_library_dev']:.2e}, rounds {r['hold_rounds']}, SD ln "
                     f"{r['sd_ln_factor']:.4f}: {'pass' if r['passes'] else 'not built'}")
    L.append("G0e: depth ratios " + ", ".join(f"{s} {r:.3f}" for s, r in res["G0e"]["depth_ratios"].items()))
    for name in ("G0e", "G0f"):
        for j, v in res[name]["axes"].items():
            if "directions" not in v:
                L.append(f"{name} {j}: {v['why']}: held")
                continue
            L.append(f"{name} {j}: " + ", ".join(f"{s} rho {d['spearman']:+.3f} (r {d['pearson']:+.3f})"
                                                 for s, d in v["directions"].items())
                     + f": {'pass' if v['passes'] else 'held'} (bar 0.3)")
    g = res["G0g"]
    L.append("G0g: " + ", ".join(f"s({s}) {x:.3f}" for s, x in g["positions"].items())
             + f": {'pass' if g['passes'] else 'S held'}")
    return "\n".join(L) + "\n"


def cmd_tables(a) -> int:
    try:
        res = make_tables(Path(a.data), Path(a.out), a.threads)
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return EXIT_REFUSED
    print(f"wrote {res['manifest']}: sha256 {res['sha256']}", flush=True)
    return 0


def cmd_g0(a) -> int:
    out = Path(a.out)
    st = stamp()
    if not st["tree_clean"] or not st["commit"]:
        print("refused: the tracked tree is not clean (or git cannot tell)", file=sys.stderr)
        return EXIT_REFUSED
    if (out / "g0.json").exists() or (out / "g0.txt").exists():
        print(f"refused: {out} holds a G0 reading already (never overwritten)", file=sys.stderr)
        return EXIT_REFUSED
    try:
        with contextlib.redirect_stdout(io.StringIO()):  # the loading code's own lines are dropped
            res = run_g0(Path(a.data), Path(a.run_dir), Path(a.decision), Path(a.manifest), a.manifest_sha256,
                         a.board, a.threads)
        txt = g0_text(res)
    except (Refused, sc.CubeRefused) as e:
        print(f"refused: {e}", file=sys.stderr)
        return EXIT_REFUSED
    except Exception:  # noqa: BLE001
        print(masked(traceback.format_exc()), file=sys.stderr, end="")
        return EXIT_CRASH
    now = stamp()
    res["provenance"] = st | {"end_commit": now["commit"], "end_clean": now["tree_clean"]}
    res["script_sha256"] = sha256(SCRIPT_PATH)
    if now["commit"] != st["commit"] or not now["tree_clean"]:
        print("refused: the tree changed while G0 ran: no reading written", file=sys.stderr)
        return EXIT_REFUSED
    out.mkdir(parents=True, exist_ok=True)
    (out / "g0.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "g0.txt").write_text(txt)
    print(txt, end="")
    return 0 if res["outcome"] == "CONTINUE" else EXIT_STOP


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("tables")
    t.add_argument("--data", required=True, help="the released files' folder (data/raw)")
    t.add_argument("--out", default=str(DEFAULT_OUT))
    t.add_argument("--threads", type=int, default=8)
    g = sub.add_parser("g0")
    g.add_argument("--data", required=True)
    g.add_argument("--run-dir", required=True, help="the frozen gata4 run folder (factors/), read only")
    g.add_argument("--decision", default=str(HERE / "results" / "gata4" / "gata4_decisions.json"))
    g.add_argument("--manifest", default=str(DEFAULT_OUT / MANIFEST))
    g.add_argument("--manifest-sha256", required=True, help="the sha256 row K31 records")
    g.add_argument("--board", default="T3:bcatenin")
    g.add_argument("--out", required=True)
    g.add_argument("--threads", type=int, default=16)
    a = ap.parse_args(argv)
    return {"tables": cmd_tables, "g0": cmd_g0}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
