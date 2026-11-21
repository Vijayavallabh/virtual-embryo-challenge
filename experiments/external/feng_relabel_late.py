"""Feng et al. 2022 (register row X002): the frozen permitted-only labeller rerun on the post-window daily stages.

Ledger row T67 (experiments/t1_m1/timing_read.py) measures each gene's late slope inside each recomputed Feng group
over E14.5, E15.5, E16.5 and E17.5. Its E14.5 point is the existing recomputed statistics (row N54's run of
experiments/external/feng_relabel.py); this script writes the other three stages' statistics, labelled by the same
labeller. Nothing of feng_relabel.py changes: this script imports its reader, its labelling (`feng_relabel.relabel`,
with the frozen yaml's every value) and the clean reduction's statistics layout (`feng2022_geo.group_stats`, the same
seeded halves). The one difference is the list of stages handed to the labeller: LATE_STAGES in place of the yaml's
`stages`. So the control genes of each marker score are drawn from the late stages' own mean expression, by the yaml's
rule, and the yaml's unassigned cap is read on the late stages.

What it reads, and nothing else: the clean reduction's cd1_legal_counts.npz, cd1_legal_genes.txt and, from
cd1_legal_obs.csv, the cell names and the columns stage, zone and strain (feng_relabel.read_obs; the authors' cell
type column is never read). Only cells whose stage is one of LATE_STAGES are labelled; no in-window stage is in the
clean reduction (E10.5-E13.5 were dropped before it was written) and E9.5 and E14.5 are not taken here.

It runs once, under its own target-free row (an N row of experiments/harness/TRIES.md naming this script), and only
in order (`preflight`, all checked before any Feng file is opened): the yaml is feng_relabel.yaml as committed at HEAD
with sha256 YAML_SHA256, frozen, no gene UNSOURCED and no choice PENDING; row N54's relabelling record is committed at
HEAD from that yaml with route recomputed; the tree is clean; X002's register row names this use (X002_TOKEN); the row
is in the ledger and names this script; and no late relabelling is recorded or started (neither its folder under data/
nor RECORD, in the tree or at HEAD). It writes data/derived/feng_relabel_late/ (cd1_stats_relabel_late.npz, its
cells table and MANIFEST.json) and a copy of the manifest to RECORD, which the main session commits at once: row T67's
read refuses until RECORD is committed and equal to the manifest, and until the row names the statistics' sha256.
A run that stops after preflight writes a stop record (usable false, no output listed) and exits STOP_EXIT; it is not
rerun under the same row. It reads no T1 cell, no truth and no score.

Run (DGX campaign checkout ~/v/vec-campaign-d, repo root, minutes of CPU):
  uv run --frozen --no-sync python -I -B experiments/external/feng_relabel_late.py --row N<id>
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng2022_geo  # noqa: E402
import feng_relabel as fr  # noqa: E402
import filterio as fio  # noqa: E402

SOURCE = "feng_relabel_late"
LATE_STAGES = ("E15.5", "E16.5", "E17.5")
YAML_SHA256 = "883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49"  # row T62's frozen settings
LATE_DIR = "data/derived/feng_relabel_late"
STATS = "cd1_stats_relabel_late.npz"
CELLS = "cd1_relabel_late_cells.csv"
RECORD = fio.HERE / "results" / "feng_relabel_late" / "MANIFEST.json"
X002_TOKEN = "feng_relabel_late"  # X002's register row must name this use (the folder) before any run
SCRIPT_NAME = "feng_relabel_late.py"
N_ROW = re.compile(r"N[0-9]+[a-z]?")
STOP_EXIT = 3
CODE = [*fr.CODE, Path(__file__)]


class Refused(RuntimeError):
    """A rule of this script stops the run before any Feng file is opened."""


def x002_names_this_use(path: Path | None = None) -> bool:
    """EXTERNAL_DATA.md's X002 row (its own line) names X002_TOKEN."""
    text = Path(fio.REPO / "EXTERNAL_DATA.md" if path is None else path).read_text(encoding="utf-8")
    return any(line.startswith("| X002 |") and X002_TOKEN in line for line in text.splitlines())


def row_line(row: str, path: Path | None = None) -> str | None:
    """The ledger line of `row` (its first cell's first word), or None."""
    for ln in Path(path or fr.TRIES).read_text(encoding="utf-8").splitlines():
        if ln[:2] == "| " and ln[2:].split("|", 1)[0].split()[:1] == [row]:
            return ln
    return None


def late_config(cfg: dict) -> dict:
    """The frozen yaml's values with LATE_STAGES as the stages handed to the labeller; nothing else differs."""
    out = copy.deepcopy(dict(cfg))
    out["stages"] = list(LATE_STAGES)
    return out


def preflight(row: str, dst: Path) -> dict:
    """Everything the run needs, checked before any Feng file is opened. Returns the gate the manifest records."""
    from vec import provenance

    if not N_ROW.fullmatch(str(row)):
        raise Refused(f"{row!r} is not an N row: the late labelling is a target-free row of its own")
    raw = fr.CONFIG.read_bytes()
    if fio.sha256_file(fr.CONFIG) != YAML_SHA256:
        raise Refused(f"feng_relabel.yaml's sha256 is not the frozen {YAML_SHA256}")
    cfg = fr.load_config(fr.CONFIG)
    if cfg.get("status") != "frozen" or fr.unsourced_genes(cfg) or fr.pending_choices(cfg):
        raise Refused("the yaml is not frozen, or holds UNSOURCED genes or PENDING choices")
    if fr.committed_at_head(fr.CONFIG) != raw:
        raise Refused("the yaml differs from its committed version at HEAD")
    rec = fr.committed_at_head(fr.RECORD)
    first = json.loads(rec) if rec is not None else {}
    if (first.get("config") or {}).get("sha256") != YAML_SHA256 or first.get("route") != "recomputed":
        raise Refused("row N54's relabelling record is not at HEAD from the frozen yaml with route recomputed")
    if not x002_names_this_use():
        raise Refused(f"EXTERNAL_DATA.md's X002 row does not name {X002_TOKEN!r}: the use is logged before the run")
    line = row_line(row)
    if line is None or SCRIPT_NAME not in line:
        raise Refused(f"{row} has no row in experiments/harness/TRIES.md naming {SCRIPT_NAME}: the row comes first")
    commit, unclean = provenance.git_state(fio.REPO)
    if unclean is None or unclean:
        raise Refused(f"the tree is not clean ({unclean}): the late labelling runs once, on a committed tree")
    if Path(dst).exists() or RECORD.exists() or fr.committed_at_head(RECORD) is not None:
        raise Refused("a late relabelling is already recorded or was started: a rerun is a new ledger row")
    return {"row": row, "commit": commit, "clean": True, "config_sha256": YAML_SHA256,
            "first_record_route": first.get("route")}


def run(cfg: dict, src: Path, dst: Path, gate: dict) -> dict:
    """Read the clean reduction's late cells (manifest-checked), label them with the frozen labeller, write the
    statistics, the cells table and the manifest. Returns the manifest."""
    if not gate or gate.get("clean") is not True:
        raise Refused("no clean gate: the late labelling runs only through main(), after preflight")
    inp = cfg["inputs"]
    paths = {k: src / inp[k] for k in ("counts", "genes", "obs")}
    entries = {k: fio.verify_manifest_output(p) for k, p in paths.items()}
    obs = fr.read_obs(paths["obs"], cfg["obs_columns"])
    if not np.all(obs["strain"] == cfg["expect_strain"]):
        raise Refused(f"a cell's strain is not {cfg['expect_strain']}")
    genes = paths["genes"].read_text().splitlines()
    C = sparse.load_npz(paths["counts"]).tocsr()
    if C.shape != (len(obs["cell"]), len(genes)):
        raise Refused(f"counts {C.shape} do not match the obs table and genes")
    lcfg = late_config(cfg)
    keep = np.flatnonzero(np.isin(obs["stage"], list(LATE_STAGES)))
    names, stage, zone = obs["cell"][keep], obs["stage"][keep], obs["zone"][keep]
    C = C[keep]
    res = fr.relabel(C, genes, stage, lcfg)
    st = feng2022_geo.group_stats(C, stage, zone, res["labels"], names, genes)
    dst.mkdir(parents=True, exist_ok=True)
    fio.save_npz_atomic(dst / STATS, **st)
    fr._write_cells(dst / CELLS, names, stage, zone, res["labels"], res["scores"])
    usable = bool(res["cap"]["passed"])
    outputs = {STATS: {"sha256": fio.sha256_file(dst / STATS), "usable": usable, "groups": int(len(st["keys"])),
                       "half_seed": feng2022_geo.HALF_SEED},
               CELLS: {"sha256": fio.sha256_file(dst / CELLS)}}
    manifest = {"source": SOURCE, "method": cfg["method"], "register_row": "X002",
                "config": {"file": fr._rel(fr.CONFIG), "sha256": fio.sha256_file(fr.CONFIG),
                           "status": cfg.get("status")},
                "stages_handed_to_the_labeller": list(LATE_STAGES), "stages": list(LATE_STAGES),
                "inputs": {k: {"file": fr._rel(p), "sha256": entries[k]["sha256"]} for k, p in paths.items()},
                "obs_columns_read": ["<cell names>", *cfg["obs_columns"]], "authors_labels_read": False,
                "unassigned_share": res["unassigned_share"], "cap": res["cap"], "usable": usable,
                "group_cells": res["group_cells"], "label_counts": res["label_counts"], "sets": res["sets"],
                "unassigned_kinds": res["unassigned_kinds"], "gate": dict(gate), "outputs": outputs,
                "provenance": fio.provenance(CODE)}
    fio.write_manifest(dst, manifest)
    return manifest


def write_record(manifest_path: Path, record: Path | None = None) -> Path:
    record = RECORD if record is None else record
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_bytes(Path(manifest_path).read_bytes())
    return record


def write_stop(dst: Path, cfg: dict, gate: dict, reason: str) -> dict:
    """A run that stopped after preflight: a manifest listing no output (no reader can take a file from it), usable
    false, and its copy at RECORD."""
    manifest = {"source": SOURCE, "method": cfg.get("method"), "register_row": "X002",
                "config": {"file": fr._rel(fr.CONFIG), "sha256": fio.sha256_file(fr.CONFIG)},
                "stages": list(LATE_STAGES), "gate": dict(gate), "stopped": True, "reason": str(reason)[:500],
                "usable": False, "authors_labels_read": False, "outputs": {}, "provenance": fio.provenance(CODE)}
    Path(dst).mkdir(parents=True, exist_ok=True)
    fio.write_manifest(dst, manifest)
    write_record(Path(dst) / "MANIFEST.json")
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--row", required=True, help="the target-free N row of experiments/harness/TRIES.md for this run")
    a = ap.parse_args(argv)
    try:  # nothing of Feng's is opened in this block
        cfg = fr.load_config(fr.CONFIG)
        src = (fio.REPO / cfg["inputs"]["dir"]).resolve()
        dst = fr.output_dir_ok(fio.REPO / LATE_DIR)
        gate = preflight(a.row, dst)
    except (Refused, fr.Refused) as e:
        raise SystemExit(f"refused: {e}") from None
    dst.mkdir(parents=True)
    (dst / fr.STARTED).write_text(json.dumps({"gate": gate}, indent=1) + "\n")
    try:
        man = run(cfg, src, dst, gate)
    except BaseException as e:  # any stop after preflight is recorded; nothing is rerun under the row
        write_stop(dst, cfg, gate, f"{type(e).__name__}: {e}")
        print(f"late labelling stopped after preflight; commit {RECORD.relative_to(fio.REPO)} now")
        if not isinstance(e, Exception):
            raise
        return STOP_EXIT
    rec = write_record(dst / "MANIFEST.json")
    print(f"late labelling written; usable {man['usable']}; commit {rec.relative_to(fio.REPO)} now")
    return 0


if __name__ == "__main__":
    sys.exit(main())
