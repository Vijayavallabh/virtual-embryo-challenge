"""Row N54 of experiments/harness/TRIES.md (T1; a diagnostic, not a try; target-free): what the permitted-only
relabelling of Feng's cells (experiments/external/feng_relabel.py, register row X002) does to M1's direction.

It runs only after the captain's ruling on the second T1 file and after the recomputed statistics file is written and
hashed: it refuses unless (1) row N54 is in TRIES.md, (2) the decision-log entry named by --ruling is in docs/PLAN.md
and names N54 anywhere in the entry, and is the ruling the relabelling's gate recorded, (3) the sha256 given by
--recomputed-sha256 is the recomputed file's and its manifest's, (4) the relabelling's manifest is committed at HEAD
as experiments/external/results/feng_relabel/MANIFEST.json, byte for byte, (5) the file passes feng_dirs'
checks (gate, frozen yaml as it stands; an unusable file is described but forms no direction), (6) the tree is clean
and (7) no output exists (a rerun is a new row). A relabelling that stopped after preflight (its stop record, route
none by terminal.on_stop) has nothing to list, and the listing says so and exits.

Per direction group (CM, Endocardium, Endothelium, Epicardium) it lists: the recomputed cells at each end of the
bracket (E9.5, E14.5); the unassigned share per stage and the cap's verdict (from the relabelling's manifest);
whether the group fell back and to what; the terminal rule's route and the one builder command it fixes; and, over
"the genes M1 moves" (a value in the present direction, and on
neither the pre-registered held-flat list nor R2's union), the Pearson correlation and the cosine between the
recomputed and the present direction, their norms, the genes whose factor log f = 0.4 x d changes sign (both versions
non-zero), the same among genes with abs(log f) >= log(1.1) in either version (external_shift's n_up/n_down bar), and
per version the counts above that bar. Beside it, deciding nothing: the recomputed groups' cells at each end under the
yaml's sensitivity settings, from the relabelling's per-cell scores (written at full precision, so a cell on a
threshold falls as in the run).

It decides nothing: nothing reads its output, and no rule here moves a file (the route is the labeller's, fixed in the
yaml before the run). It never reads the authors' label column, but it does use the authors' labels once, indirectly:
the present direction it compares against is formed on the authors' groups (the clean reduction's cd1_stats.npz,
through feng_dirs, as every build reads it). The per-group agreement with the authors' labels is a separate script
(n54_authors_agreement.py), run once afterwards for the method report only. It reads Feng's two statistics files
(manifest-checked), the relabelling's manifest and per-cell scores, the gene names of the released
T1/E10.5_RNA.h5ad (names only, the gene order every T1 build uses), the pre-registered held-flat list and R2's
pinned union; no target, score, board, built file or T1 cell count.

Writes experiments/t1_m1/results/n54_relabel_listing/n54.{json,txt}.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in ("t1_hygiene", "t1_direction", "external"):
    sys.path.insert(0, str(REPO / "experiments" / _p))
sys.path.insert(0, str(HERE))
import feng_dirs  # noqa: E402
import feng_relabel  # noqa: E402
import filterio as fio  # noqa: E402

ruling_logged = feng_relabel.ruling_logged  # the decision-log entry, read whole (feng_relabel's gate uses the same)

ROW = "N54"
GROUPS = ("CM", "Endocardium", "Endothelium", "Epicardium")  # vec.t1_state_transport.R1_FENG_GROUPS
FRACTION = 0.4  # feng_dirs.bracket_fraction(10.5, 12.5), the step M1 applies
BAR = float(np.log(1.1))  # vec.t1_scope.external_shift's n_up / n_down threshold
T1_FILE = REPO / "data" / "raw" / "T1" / "E10.5_RNA.h5ad"
OUT = HERE / "results" / "n54_relabel_listing"
PLAN = REPO / "docs" / "PLAN.md"


# ------------------------------------------------------------------------------------------------ pure parts


def moved_mask(d_present: np.ndarray, genes: Iterable[str], held: Iterable[str]) -> np.ndarray:
    """The genes M1 moves in the present direction: a finite value, and on no held list."""
    hs = set(map(str, held))
    return np.isfinite(np.asarray(d_present, dtype=np.float64)) & ~np.array([str(g) in hs for g in genes], dtype=bool)


def compare(d_new: np.ndarray, d_old: np.ndarray, moved: np.ndarray, fraction: float = FRACTION,
            bar: float = BAR) -> dict:
    """The recomputed direction against the present one over the moved genes. A gene the recomputed direction has no
    value for counts as 0 there (external_shift's nan_to_num), as it would in a build."""
    old = np.asarray(d_old, dtype=np.float64)[moved]
    new = np.nan_to_num(np.asarray(d_new, dtype=np.float64)[moved], nan=0.0)
    f_old, f_new = fraction * old, fraction * new
    out = {"genes_moved": int(moved.sum()),
           "genes_without_recomputed_value": int((~np.isfinite(np.asarray(d_new, dtype=np.float64)[moved])).sum())}
    if out["genes_moved"] < 2:
        return {**out, "pearson": None, "cosine": None}
    n_old, n_new = float(np.linalg.norm(old)), float(np.linalg.norm(new))
    out["norm_present"], out["norm_recomputed"] = n_old, n_new
    out["cosine"] = float(old @ new / (n_old * n_new)) if n_old > 0 and n_new > 0 else None
    out["pearson"] = float(np.corrcoef(old, new)[0, 1]) if np.std(old) > 0 and np.std(new) > 0 else None
    flip = (np.sign(f_old) * np.sign(f_new)) < 0
    large = (np.abs(f_old) >= bar) | (np.abs(f_new) >= bar)
    out["sign_changes"] = int(flip.sum())
    out["sign_changes_large"] = int((flip & large).sum())
    out["n_up_present"], out["n_down_present"] = int((f_old > bar).sum()), int((f_old < -bar).sum())
    out["n_up_recomputed"], out["n_down_recomputed"] = int((f_new > bar).sum()), int((f_new < -bar).sum())
    return out


def group_rows(manifest: Mapping, record: Mapping, groups: Iterable[str] = GROUPS) -> dict:
    """Per group, from the relabelling's manifest and the recomputed direction's record: the cells at each end, whether
    a fallback fired and what the group took. Shares per stage are the manifest's own."""
    fb = record.get("fallbacks", {})
    return {g: {"cells": dict(manifest["group_cells"].get(g, {})),
                "fallback_fired": bool(fb.get(g, {}).get("fired", False)),
                "takes": fb.get(g, {}).get("takes")} for g in groups}


def sensitivity(scores: Mapping[str, np.ndarray], stage: np.ndarray, cfg: Mapping) -> list[dict]:
    """The recomputed groups' cells at each end under each of the yaml's sensitivity settings (deciding nothing)."""
    stages = list(cfg["stages"])
    rows = []
    for s in cfg.get("sensitivity", []):
        labels = feng_relabel.assign(scores, cfg, s["min_score"], s["min_margin"])
        rows.append({"setting": dict(s), "group_cells": feng_relabel.group_cells(labels, stage, stages, cfg["groups"]),
                     "unassigned_share": feng_relabel.unassigned_shares(labels, stage, stages,
                                                                         cfg["assign"]["unassigned"])})
    return rows


def read_scores(path: Path, cfg: Mapping) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """The relabelling's per-cell scores and stages (its cells table; no authors' column exists there)."""
    with open(path, newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        rows = list(rd)
    i_stage = header.index("stage")
    stage = np.array([r[i_stage] for r in rows], dtype=str)
    scores = {n: np.array([float(r[header.index(f"score:{n}")]) for r in rows]) for n in cfg["sets"]}
    return scores, stage


def stopped(data_manifest: Path) -> str | None:
    """The reason to list nothing: the relabelling stopped after preflight (feng_relabel's stop record; its route is
    none by terminal.on_stop and it wrote no statistics to list), or it has not run. None when there is a run to
    list."""
    p = Path(data_manifest)
    if not p.exists():
        return "the permitted-only relabelling has not written its manifest: there is nothing to list"
    man = json.loads(p.read_text())
    if man.get("stopped"):
        return (f"the relabelling stopped after preflight ({man.get('reason')}); its route is {man.get('route')} by "
                "terminal.on_stop, and N54 lists nothing")
    return None


def record_committed(data_manifest: Path, record: Path | None = None) -> bool:
    """The relabelling's manifest is committed at HEAD as its record, byte for byte."""
    head = feng_relabel.committed_at_head(feng_relabel.RECORD if record is None else record)
    return head is not None and head == Path(data_manifest).read_bytes()


def text(rep: Mapping) -> str:
    lines = [rep["label"], f"cap: {rep['cap']}", f"route: {rep['route']} ({' '.join(rep['command'] or [])})"]
    for g, r in rep["groups"].items():
        c = r["compare"]
        lines.append(f"{g}: cells {r['cells']}; fallback {r['fallback_fired']} -> {r['takes']}; moved genes "
                     f"{c['genes_moved']}; pearson {c.get('pearson')}; cosine {c.get('cosine')}; sign changes "
                     f"{c.get('sign_changes')} (large {c.get('sign_changes_large')})")
    lines.append("Decides nothing. The authors' label column is not read here; the present direction compared against "
                 "is formed on the authors' groups.")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ run


def main(argv=None) -> None:
    import hygiene
    import n12_hygiene_size as n12

    from vec import provenance, t1_scope

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ruling", required=True, help="the docs/PLAN.md decision-log entry of the captain's ruling")
    ap.add_argument("--recomputed-sha256", required=True, help="the recomputed statistics file's sha256, as hashed")
    ap.add_argument("--t1-file", default=str(T1_FILE))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if ROW not in n12.tries_rows():
        sys.exit(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if not ruling_logged(a.ruling, ROW, PLAN):
        sys.exit(f"docs/PLAN.md has no decision-log entry '{a.ruling}' naming {ROW}: the listing follows the ruling")
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{ROW} has already written its outputs; a rerun is a new row (experiments/harness/TRIES.md)")
    start = provenance.git_state(provenance.REPO_ROOT)
    if start[1] is None or start[1]:
        sys.exit(f"the tree is not clean ({start[1]}): {ROW} runs once, on a committed tree")
    rcfg = feng_dirs.load_relabel_config()
    rpath = feng_dirs.recomputed_stats_path(rcfg)
    stop = stopped(rpath.parent / "MANIFEST.json")
    if stop is not None:
        sys.exit(stop)
    entry = fio.verify_manifest_output(rpath)
    if entry["sha256"] != a.recomputed_sha256:
        sys.exit("--recomputed-sha256 is not the recomputed file's sha256 in its manifest")
    man_path = rpath.parent / "MANIFEST.json"
    if not record_committed(man_path):
        sys.exit(f"the relabelling's manifest is not committed as {feng_relabel.RECORD.relative_to(REPO)}: commit it "
                 "first, so the run is on the record before anything reads it")
    try:  # the gate, the frozen yaml as it stands; an unusable file is described but forms no direction
        _, manifest = feng_dirs.check_recomputed(rpath, require_usable=False)
    except ValueError as e:
        sys.exit(f"refused: {e}")
    if manifest["gate"]["ruling"] != a.ruling:
        sys.exit(f"--ruling {a.ruling!r} is not the ruling the relabelling ran under ({manifest['gate']['ruling']!r})")
    genes = n12.t1_gene_names(Path(a.t1_file))
    present, rec_present = feng_dirs.feng_directions(genes)
    usable = bool(entry.get("usable"))
    if usable:
        recomputed, rec_new = feng_dirs.feng_directions(genes, groups="recomputed")
    else:  # the cap failed: no direction is formed; the listing still reports the cells and shares
        recomputed, rec_new = {}, {"fallbacks": {}}
    held = set(t1_scope.held_flat_genes()) | set(hygiene.load_union()["genes"])
    groups = {}
    for g, row in group_rows(manifest, rec_new).items():
        if g in present and g in recomputed:
            row["compare"] = compare(recomputed[g], present[g], moved_mask(present[g], genes, held))
        else:
            row["compare"] = {"genes_moved": 0, "present": g in present, "recomputed": g in recomputed}
        groups[g] = row
    scores, stage = read_scores(rpath.parent / rcfg["output"]["cells"], rcfg)
    route = manifest["route"]
    report = {"label": "N54: the permitted-only relabelling of Feng's cells against M1's present direction, which is "
                       "formed on the authors' groups (target-free; decides nothing; the authors' label column "
                       "not read)",
              "ruling": a.ruling, "usable": usable, "cap": manifest["cap"],
              "route": route, "command": feng_relabel.build_command(route, rcfg),
              "unassigned_share": manifest["unassigned_share"], "unassigned_kinds": manifest["unassigned_kinds"],
              "groups": groups,
              "sensitivity": sensitivity(scores, stage, rcfg), "fraction": FRACTION, "bar": BAR,
              "inputs": {"recomputed": {"file": str(rpath.relative_to(REPO)), "sha256": entry["sha256"]},
                         "relabel_config_sha256": manifest["config"]["sha256"],
                         "present": {"file": rec_present["file"], "sha256": rec_present["sha256"]},
                         "t1_file": str(Path(a.t1_file).resolve().relative_to(REPO)),
                         "t1_file_sha256": hygiene.sha256_file(a.t1_file)},
              "stamp": {"commit": start[0], "tree_clean_at_start": True, "tree_clean_at_end": None,
                        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    end_commit, end_unclean = provenance.git_state(provenance.REPO_ROOT)
    report["stamp"]["tree_clean_at_end"] = None if end_unclean is None else (not end_unclean and end_commit == start[0])
    out.mkdir(parents=True, exist_ok=True)
    (out / "n54.json").write_text(json.dumps(report, indent=1, default=float) + "\n")
    (out / "n54.txt").write_text(text(report))
    print(text(report), end="")


if __name__ == "__main__":
    main()
