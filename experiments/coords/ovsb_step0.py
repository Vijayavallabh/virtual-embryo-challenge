#!/usr/bin/env python
"""Step 0 of the file T2embryo_s2x_ovsb.h5ad (row R48 of experiments/harness/TRIES.md, the captain's order of
2026-11-24; this step is row E56, one embryo try by its check (ii)): three checks of the builder's
`--expression series_back` on the released brackets E7.5 and E8.0 alone. They gate the file: any STOP means no file, and
the handed-over files stand. No number here comes from a test score, and nothing here reads a test board's counts, a
test file or a target stage.

The checks (STATEMENT, fixed before any data were read; written into the reading, never tuned):
  (i)   Reported only: the factor's log SD, its count of genes beyond twofold, the matched labels' share of E8.0's
        cells.
  (ii)  The carry check: apply the factor to E8.0's cells of the embryonic labels that both brackets hold with at least
        100 cells; STOP if it raises the cell-weighted SD over genes of log(E8.0 mean proportions / E7.5's).
  (iii) STOP unless the rendered E8.0 cells' median depth lies within 10% of E7.5's.

How each is computed (METHOD, fixed with this script before it ever ran on data):
  The brackets are loaded as the builder loads them: data/raw/T2/E7.5.h5ad and E8.0.h5ad through vec.prep.load_stage
  with the embryo test board's gene panel (the public gene list, read by code) and drop_duplicates=True. The factor and
  its matched extraembryonic labels are vec.t2_build.series_back_factor's, the ones the build applies (SERIES_RULE's
  gene_factor from E8.0 to E7.5, weights E8.0's label shares).
  (i)   the SD over the panel genes of log(factor) (population SD, as the build record's sd_log_factor); the genes with
        |log factor| above log 2; the share of E8.0's cells whose label is a matched label.
  (ii)  The embryonic labels are the labels both brackets hold with at least SERIES_MIN_CELLS (100) cells each that are
        not matched extraembryonic labels. For each such label L: p75 and p80, the mean within-cell panel proportions
        of L's cells at E7.5 and at E8.0 (series_factor's own proportions), and p80f, the same at E8.0 with the factor
        applied; s(L) = the SD over the panel genes of log((p80 + eps) / (p75 + eps)), eps = SERIES_EPS (the factor's
        own). The figure is the mean of s(L) over the labels weighted by E8.0's cells of each label, without and with
        the factor. The factor is applied in expectation, as the build's count-space step applies it (vec.core
        apply_gene_factor): a cell's panel counts c become phi c where phi < 1 (binomial thinning's mean), c + P n
        (phi - 1) where phi > 1 (the Poisson top-up's mean; P the cell's kNN-smoothed panel proportions over the E8.0
        bracket, RENDER_K neighbours, as the build smooths them; n its panel counts) and c where phi = 1, and its
        proportions are those counts over their sum. In expectation, so that the check carries no sampling noise: the
        build's own draws would add thinning and top-up noise to every label's mean proportions, which raises the SD
        whatever the factor does. The thinning to E7.5's depth that follows changes no expected proportion and is left
        out. STOP (ii) if the figure with the factor is above the figure without it; STOP (ii) as well when no
        embryonic label qualifies (the check cannot pass).
  (iii) Every E8.0 cell rendered by vec.t2_build.series_back_expression exactly as the build renders its drawn E8.0
        rows (all of E8.0's cells in place of the draw, the build's seed and stream [BUILD_SEED_OFFSET + 1, 7]); the
        rendered cells' median counts per cell (their library sizes) against E7.5's median counts per cell (over all
        its cells, as series_back reads it). STOP (iii) unless |rendered - E7.5's| <= 0.10 x E7.5's.
  The outcome line is exactly "Step 0: PASS", "Step 0: STOP (ii)", "Step 0: STOP (iii)" or "Step 0: STOP (ii, iii)".

Refused: unless row E56 is in experiments/harness/TRIES.md and the tree is clean (vec.provenance's rule: its row is
committed before the step runs); when the results folder already holds a reading (a reading is never overwritten);
when vec is not loaded from this checkout's src/. It writes reading.txt (the outcome line, then the figures, the
statement, the method and the code that made it) and reading.json into the results folder it is given, and prints only
"wrote <path>". The launcher experiments/builds/embryo_ovsb.sh builds the file only on a committed reading of this
script at experiments/coords/results/ovsb_step0/ whose outcome is "Step 0: PASS".

Usage (repo root, after row E56 is committed; BLAS threads capped as docs/COMPUTE.md says):
  uv run --frozen --no-sync python -I -B experiments/coords/ovsb_step0.py --out experiments/coords/results/ovsb_step0
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

CODE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CODE_ROOT / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
from scipy import sparse  # noqa: E402

from vec import harness  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.prep import load_stage  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

REPO = CODE_ROOT
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
ROW = "E56"
RULING = "R48"
FILE = "T2embryo_s2x_ovsb.h5ad"
BOARD = "T2:embryo:test_interp"  # its gene panel only
STAGES = {"before": "E7.5.h5ad", "after": "E8.0.h5ad"}
RESULTS = "experiments/coords/results/ovsb_step0"  # where the launcher reads the committed reading
SEED = [harness.BUILD_SEED_OFFSET + 1, 7]  # the build's seed and expression stream
DEPTH_BOUND = 0.10
OUTCOMES = ("Step 0: PASS", "Step 0: STOP (ii)", "Step 0: STOP (iii)", "Step 0: STOP (ii, iii)")
STATEMENT = ("Step 0 of T2embryo_s2x_ovsb.h5ad (the file of row R48; this step is row E56), fixed before any data were "
             "read and not to be tuned. (i) Reported only: the factor's log SD, its count of genes beyond twofold, "
             "the matched labels' share of E8.0's cells. (ii) The carry check: apply the factor to E8.0's cells of the "
             "embryonic labels that both brackets hold with at least 100 cells; STOP if it raises the cell-weighted "
             "SD over genes of log(E8.0 mean proportions / E7.5's). (iii) STOP unless the rendered E8.0 cells' median "
             "depth lies within 10% of E7.5's. Any STOP: no file, and the handed-over files stand.")
METHOD = __doc__[__doc__.index("How each is computed"):__doc__.index("  The outcome line")].strip()


# ------------------------------------------------------------------------------------------------ the checks
def mean_props(C, rows) -> np.ndarray:
    """The mean within-cell panel proportions of `rows` (series_factor's own proportions)."""
    sub = sparse.csr_matrix(C[rows], dtype=np.float64)
    tot = np.asarray(sub.sum(1)).ravel()
    return np.asarray((sparse.diags(1.0 / np.maximum(tot, 1.0)) @ sub).mean(0)).ravel()


def expected_counts(C, phi, P) -> np.ndarray:
    """vec.core.apply_gene_factor's expected counts: phi c where phi < 1, c + P n (phi - 1) where phi > 1 (n the cell's
    panel counts), c where phi = 1. C cells x genes, P the cells' smoothed proportions (rows summing to 1)."""
    C = np.asarray(C, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    out = C.copy()
    down, up = phi < 1.0, phi > 1.0
    out[:, down] = C[:, down] * phi[down]
    out[:, up] = C[:, up] + np.asarray(P, dtype=np.float64)[:, up] * C.sum(1, keepdims=True) * (phi[up] - 1.0)
    return out


def embryonic_labels(names75, names80, matched) -> list[str]:
    """The labels both brackets hold with at least SERIES_MIN_CELLS cells each, less the matched extraembryonic ones."""
    n75, n80 = np.asarray(names75).astype(str), np.asarray(names80).astype(str)
    both = sorted(set(n75) & set(n80))
    return [lab for lab in both if lab not in set(matched) and (n75 == lab).sum() >= tb.SERIES_MIN_CELLS
            and (n80 == lab).sum() >= tb.SERIES_MIN_CELLS]


def carry_check(before, after, genes, phi, matched) -> dict:
    """(ii): the cell-weighted SD over genes of log(E8.0 mean proportions / E7.5's) over the embryonic labels, without
    and with the factor applied in expectation (METHOD). raised: with > without, or no embryonic label."""
    _, C75, _, _ = tb._stage_counts(before, genes)
    E80, C80, lib80, _ = tb._stage_counts(after, genes)
    n75 = np.asarray(before.obs["celltype"]).astype(str)
    n80 = np.asarray(after.obs["celltype"]).astype(str)
    labels = embryonic_labels(n75, n80, matched)
    if not labels:
        return {"labels": [], "raised": True, "reason": "no embryonic label held by both brackets with at least "
                f"{tb.SERIES_MIN_CELLS} cells: the carry check cannot pass"}
    rows = np.flatnonzero(np.isin(n80, labels))
    S = tb._smoothed_props(E80, C80, lib80, rows, tb.RENDER_K)  # as the build smooths a cell's profile
    P = S / np.maximum(S.sum(1, keepdims=True), 1e-12)
    C0 = np.asarray(C80[rows].toarray(), dtype=np.float64)
    Ef = expected_counts(C0, phi, P)
    # each cell's proportions without and with the factor, by the same arithmetic (a factor of 1 gives the same bits)
    P0 = C0 / np.maximum(C0.sum(1, keepdims=True), 1.0)
    Pf = Ef / np.maximum(Ef.sum(1, keepdims=True), 1e-300)
    where = {lab: np.flatnonzero(n80[rows] == lab) for lab in labels}
    eps = tb.SERIES_EPS
    per = []
    for lab in labels:
        p75 = mean_props(C75, np.flatnonzero(n75 == lab))
        p80 = P0[where[lab]].mean(0)
        p80f = Pf[where[lab]].mean(0)
        per.append({"label": lab, "weight": float(len(where[lab])),
                    "sd_without": float(np.log((p80 + eps) / (p75 + eps)).std()),
                    "sd_with": float(np.log((p80f + eps) / (p75 + eps)).std())})
    w = np.array([r["weight"] for r in per])
    w /= w.sum()
    for r, wi in zip(per, w, strict=True):
        r["weight"] = float(wi)  # E8.0's share of the label among the embryonic labels' cells (no count is kept)
    without = float(sum(wi * r["sd_without"] for wi, r in zip(w, per, strict=True)))
    with_ = float(sum(wi * r["sd_with"] for wi, r in zip(w, per, strict=True)))
    return {"labels": labels, "per_label": per, "sd_without": without, "sd_with": with_, "raised": with_ > without}


def depth_within(rendered: float, d75: float, bound: float = DEPTH_BOUND) -> bool:
    """(iii): the rendered median within `bound` x E7.5's of E7.5's, the bound included."""
    return bool(abs(rendered - d75) <= bound * d75)


def outcome(stop_ii: bool, stop_iii: bool) -> str:
    stops = [k for k, bad in (("ii", stop_ii), ("iii", stop_iii)) if bad]
    return "Step 0: PASS" if not stops else f"Step 0: STOP ({', '.join(stops)})"


def reading(before, after, genes) -> dict:
    """The three checks on the two brackets (loaded as the builder loads them)."""
    phi, matched = tb.series_back_factor(before, after, genes)
    n80 = np.asarray(after.obs["celltype"]).astype(str)
    log_phi = np.log(phi)
    i = {"sd_log_factor": float(log_phi.std()), "genes_beyond_2x": int((np.abs(log_phi) > np.log(2.0)).sum()),
         "n_genes": len(genes), "matched_labels": matched, "matched_share_e80": float(np.isin(n80, matched).mean())}
    ii = carry_check(before, after, genes, phi, matched)
    every = SimpleNamespace(rows=np.arange(after.n_obs), stage=np.ones(after.n_obs, dtype=np.int64))
    _, rec = tb.series_back_expression(before, after, every, genes, 0.5, SEED)
    if rec["labels"] != matched or rec["sd_log_factor"] != i["sd_log_factor"]:  # one factor, the build's
        raise SystemExit("the rendering's factor is not the factor read here: refused, nothing written")
    d75, rendered = rec["median_depth"]["before"], rec["rendered_median"]
    iii = {"rendered_median_e80": rendered, "median_e75": d75, "median_e80": rec["median_depth"]["after"],
           "ratio_to_e75": rendered / d75, "bound": DEPTH_BOUND, "within": depth_within(rendered, d75),
           "thinning_ratio": rec["ratio"], "seed": SEED}
    return {"outcome": outcome(ii["raised"], not iii["within"]), "row": ROW, "ruling": RULING, "file": FILE,
            "statement": STATEMENT, "method": METHOD, "i": i, "ii": ii, "iii": iii}


def text(res: dict) -> str:
    i, ii, iii = res["i"], res["ii"], res["iii"]
    L = [res["outcome"],
         f"(i) reported only: the factor's log SD {i['sd_log_factor']:.4f}; genes beyond twofold "
         f"{i['genes_beyond_2x']} of {i['n_genes']}; the matched labels' share of E8.0's cells "
         f"{i['matched_share_e80']:.4f} (labels: {', '.join(i['matched_labels'])})"]
    if ii["labels"]:
        L.append(f"(ii) carry check over {len(ii['labels'])} embryonic labels: the cell-weighted SD over genes of "
                 f"log(E8.0 mean proportions / E7.5's) {ii['sd_without']:.4f} without the factor, {ii['sd_with']:.4f} "
                 f"with it: {'raised (STOP)' if ii['raised'] else 'not raised'}")
        L += [f"     {r['label']}: E8.0 share {r['weight']:.4f}, SD {r['sd_without']:.4f} without, {r['sd_with']:.4f} "
              "with" for r in ii["per_label"]]
    else:
        L.append(f"(ii) carry check: {ii['reason']} (STOP)")
    L.append(f"(iii) depth: the rendered E8.0 cells' median counts per cell {iii['rendered_median_e80']:.2f} against "
             f"E7.5's {iii['median_e75']:.2f} (ratio {iii['ratio_to_e75']:.4f}, bound {iii['bound']:.0%}): "
             f"{'within' if iii['within'] else 'outside (STOP)'}")
    L += ["", "Statement (fixed before any data were read): " + res["statement"], "", res["method"], ""]
    p = res.get("provenance", {})
    L.append(f"Code: commit {p.get('commit')}, clean {p.get('clean')}, script sha256 {p.get('script_sha256')}, "
             f"t2_build sha256 {p.get('t2_build_sha256')}")
    return "\n".join(L) + "\n"


# ------------------------------------------------------------------------------------------------ the run
def has_row(path: Path | None = None) -> bool:
    """Whether the ledger holds row E56 (a table row whose first cell starts with its ID)."""
    return any(ln.startswith(f"| {ROW} ") for ln in (path or TRIES).read_text(encoding="utf-8").splitlines())


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _own_code() -> None:
    src = (CODE_ROOT / "src").resolve()
    for mod in (harness, tb):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help=f"the results folder (the launcher reads {RESULTS})")
    ap.add_argument("--data", default=str(harness.DATA_DIR), help="the released files' root (T2/ under it)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    files = {"txt": out / "reading.txt", "json": out / "reading.json"}
    existing = [str(f) for f in files.values() if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)} (a reading is never overwritten)")
    _own_code()
    if not has_row():
        sys.exit(f"row {ROW} is not in experiments/harness/TRIES.md: its row is committed before step 0 runs")
    state = repo_state(REPO)
    if state.unclean != ():
        sys.exit(f"the tree is not clean: step 0 runs on a committed tree, with row {ROW} committed")
    genes = list(load_boards()[BOARD].genes)
    paths = {k: Path(a.data) / "T2" / f for k, f in STAGES.items()}
    stages = {k: load_stage(p, genes=genes, drop_duplicates=True) for k, p in paths.items()}
    res = reading(stages["before"], stages["after"], genes)
    res["provenance"] = {"commit": state.commit, "clean": True, "script_sha256": _sha(Path(__file__)),
                         "t2_build_sha256": _sha(Path(tb.__file__)),
                         "inputs": {k: {"file": f"T2/{STAGES[k]}", "sha256": _sha(p)} for k, p in paths.items()}}
    out.mkdir(parents=True, exist_ok=True)
    files["json"].write_text(json.dumps(res, indent=1, default=float))
    files["txt"].write_text(text(res))
    for f in files.values():
        print(f"wrote {f}", flush=True)


if __name__ == "__main__":
    main()
