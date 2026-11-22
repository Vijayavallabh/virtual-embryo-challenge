"""Feng et al. 2022 (register row X002): the atlas's permitted E18 cells labelled with row N54's frozen control sets.

T1 candidate rc18 (the outside heart atlas's direction over E9.5 to E18; ledger rows N58 and T69) needs recomputed
groups at Feng's E18 stage. Rerunning experiments/external/feng_relabel.py with E18 among its stages would change
every cell's control genes (gene_means averages over the yaml's stages, and control_sets ranks genes by those means),
so it could relabel E9.5 and E14.5 cells and move row R13's route. This script does not: it labels the E18 cells alone
with the control gene sets row N54's run drew and recorded (its committed RECORD, `sets.<set>.controls`), the marker
sets and the thresholds of the frozen yaml (sha256 YAML_SHA256), and the yaml's assignment rule
(`feng_relabel.assign`); the scores are `feng_relabel.set_score` on each cell's own log1p CP10k. Nothing of
feng_relabel.py changes, and no control gene is drawn here.

The E9.5 and E14.5 statistics are row N54's file (data/derived/feng_relabel_c/cd1_stats_relabel.npz, sha256
N54_STATS_SHA256), copied unchanged: the output statistics file holds N54's rows exactly as N54 wrote them, then the E18
rows (`feng2022_geo.group_stats` on the E18 cells, the same seeded-halves rule). Before anything is written the copied
rows are checked byte for byte against N54's file (`copied_rows_identical`), and N54's file is hashed again after the
write; any difference stops the run.

What it reads, and nothing else: N54's RECORD as committed at HEAD (control sets, the set genes found, the inputs'
sha256), N54's statistics file, and the clean reduction's cd1_legal_counts.npz, cd1_legal_genes.txt and, from
cd1_legal_obs.csv, the cell names and the columns stage, zone and strain (feng_relabel.read_obs; the authors' cell type
column is never read). The three inputs must be the files N54 read (their sha256 equal to N54's record). Only cells
whose stage is E18 are labelled. It reads no T1 cell, no truth and no score.

What it writes: data/derived/feng_relabel_e18/ (cd1_stats_relabel_e18.npz, cd1_relabel_e18_cells.csv, MANIFEST.json)
and a copy of the manifest at RECORD, which the main session commits at once. It never writes feng_relabel.RECORD or
data/derived/feng_relabel_c/.

It runs once, under its own target-free N row of experiments/harness/TRIES.md naming this script, and only in order
(`preflight`, all checked before any Feng file is opened): the yaml is feng_relabel.yaml as committed at HEAD with
sha256 YAML_SHA256, frozen, no gene UNSOURCED and no choice PENDING; row N54's RECORD is committed at HEAD from that
yaml with route recomputed, usable, its statistics sha256 N54_STATS_SHA256 and a control set for every set of the yaml;
EXTERNAL_DATA.md's X002 row names X002_TOKEN; the row is in the ledger and names this script; the tree is clean; and no
E18 labelling is recorded or started (neither its folder under data/ nor RECORD, in the tree or at HEAD). A run that
stops after preflight writes a stop record (usable false, no output listed) and exits STOP_EXIT; it is not rerun under
the same row.

Run (a DGX checkout holding data/derived/feng_relabel_c/, repo root, minutes of CPU), through the launcher
experiments/t1_m1/rc18.sh (steps label-tests and label), or directly:
  uv run --frozen --no-sync python -I -B experiments/external/feng_relabel_e18.py --row N<id>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng2022_geo  # noqa: E402
import feng_relabel as fr  # noqa: E402
import filterio as fio  # noqa: E402

SOURCE = "feng_relabel_e18"
STAGE = "E18"  # Feng's own label for its last embryonic stage, read as day 18.0 (row N20b; harness.LP_BRACKET)
STAGES = (STAGE,)
COPIED_STAGES = ("E9.5", "E14.5")  # row N54's stages, copied unchanged from its statistics file
YAML_SHA256 = "883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49"  # row N54's frozen settings
N54_STATS_SHA256 = "36d1f7edcbbf3d7a80c68b51d5624a94bb753668079fe3d69aed720c71baa7e9"  # row N54's statistics
E18_DIR = "data/derived/feng_relabel_e18"
STATS = "cd1_stats_relabel_e18.npz"
CELLS = "cd1_relabel_e18_cells.csv"
RECORD = fio.HERE / "results" / "feng_relabel_e18" / "MANIFEST.json"
X002_TOKEN = "feng_relabel_e18"  # X002's register row must name this use (the folder) before any run
SCRIPT_NAME = "feng_relabel_e18.py"
N_ROW = re.compile(r"N[0-9]+[a-z]?")
STOP_EXIT = 3
STAT_FIELDS = ("keys", "genes", "n", "S_log", "S_cp", "S_nz", "depth_sum")  # feng2022_geo.group_stats' fields
ROW_FIELDS = ("n", "S_log", "S_cp", "S_nz", "depth_sum")  # one row per key
CODE = [*fr.CODE, Path(__file__)]


class Refused(RuntimeError):
    """A rule of this script stops the run (before any Feng file is opened, or, after preflight, as a recorded stop)."""


def n54_stats_path() -> Path:
    """Row N54's statistics file, as the frozen yaml names it (output.dir / output.stats)."""
    cfg = fr.load_config(fr.CONFIG)
    return fio.REPO / cfg["output"]["dir"] / cfg["output"]["stats"]


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


def n54_record(raw: bytes | None, cfg: Mapping) -> dict:
    """Row N54's RECORD (bytes as committed at HEAD), checked: from the frozen yaml, route recomputed, usable, its
    statistics sha256 N54_STATS_SHA256, and for every set of the yaml its genes found and a non-empty control list."""
    rec = json.loads(raw) if raw else {}
    if (rec.get("config") or {}).get("sha256") != YAML_SHA256 or rec.get("route") != "recomputed":
        raise Refused("row N54's relabelling record is not at HEAD from the frozen yaml with route recomputed")
    if rec.get("usable") is not True or rec.get("stopped"):
        raise Refused("row N54's relabelling record is not usable")
    stats_name = cfg["output"]["stats"]
    if ((rec.get("outputs") or {}).get(stats_name) or {}).get("sha256") != N54_STATS_SHA256:
        raise Refused(f"row N54's record does not name its statistics sha256 {N54_STATS_SHA256}")
    sets = rec.get("sets") or {}
    for name in cfg["sets"]:
        s = sets.get(name) or {}
        if not s.get("found") or not s.get("controls"):
            raise Refused(f"row N54's record has no genes found or no control set for set {name}")
    if rec.get("stages") != list(COPIED_STAGES):
        raise Refused(f"row N54's record holds the stages {rec.get('stages')}, not {list(COPIED_STAGES)}")
    return rec


def preflight(row: str, dst: Path) -> dict:
    """Everything the run needs, checked before any Feng file is opened. Returns the gate the manifest records."""
    from vec import provenance

    if not N_ROW.fullmatch(str(row)):
        raise Refused(f"{row!r} is not an N row: the E18 labelling is a target-free row of its own")
    raw = fr.CONFIG.read_bytes()
    if fio.sha256_file(fr.CONFIG) != YAML_SHA256:
        raise Refused(f"feng_relabel.yaml's sha256 is not the frozen {YAML_SHA256}")
    cfg = fr.load_config(fr.CONFIG)
    if cfg.get("status") != "frozen" or fr.unsourced_genes(cfg) or fr.pending_choices(cfg):
        raise Refused("the yaml is not frozen, or holds UNSOURCED genes or PENDING choices")
    if fr.committed_at_head(fr.CONFIG) != raw:
        raise Refused("the yaml differs from its committed version at HEAD")
    n54_record(fr.committed_at_head(fr.RECORD), cfg)
    if not x002_names_this_use():
        raise Refused(f"EXTERNAL_DATA.md's X002 row does not name {X002_TOKEN!r}: the use is logged before the run")
    line = row_line(row)
    if line is None or SCRIPT_NAME not in line:
        raise Refused(f"{row} has no row in experiments/harness/TRIES.md naming {SCRIPT_NAME}: the row comes first")
    commit, unclean = provenance.git_state(fio.REPO)
    if unclean is None or unclean:
        raise Refused(f"the tree is not clean ({unclean}): the E18 labelling runs once, on a committed tree")
    if Path(dst).exists() or RECORD.exists() or fr.committed_at_head(RECORD) is not None:
        raise Refused("an E18 labelling is already recorded or was started: a rerun is a new ledger row")
    return {"row": row, "commit": commit, "clean": True, "config_sha256": YAML_SHA256,
            "n54_record_sha256": sha256_bytes(fr.committed_at_head(fr.RECORD))}


def sha256_bytes(raw: bytes | None) -> str | None:
    return None if raw is None else hashlib.sha256(raw).hexdigest()


# ------------------------------------------------------------------------------------------------ the labelling


def fixed_sets(genes: Sequence[str], cfg: Mapping, record_sets: Mapping) -> tuple[dict, dict, dict]:
    """Per set of the yaml: its gene columns (the yaml's genes present once in Feng's symbols, in the yaml's order) and
    its control columns (row N54's recorded control genes, by name). Refused unless every set finds at least the yaml's
    min_genes_found, finds exactly the genes N54's record lists, and every recorded control gene is present once in
    Feng's symbols. Returns (set columns, control columns, the record's found and controls per set)."""
    uidx = fr.unique_index(list(genes))
    need = int(cfg["score"]["min_genes_found"])
    set_idx, ctrl_idx, shown = {}, {}, {}
    for name, s in cfg["sets"].items():
        found = [g for g in s["genes"] if g in uidx]
        if len(found) < need:
            raise Refused(f"set {name}: {len(found)} of its genes are in Feng's symbols (unique); the yaml needs "
                          f"{need}")
        rec = record_sets.get(name) or {}
        if list(rec.get("found") or []) != found:
            raise Refused(f"set {name}: the genes found differ from row N54's record")
        controls = list(rec.get("controls") or [])
        absent = [g for g in controls if g not in uidx]
        if not controls or absent:
            raise Refused(f"set {name}: row N54's control genes are missing or not unique in Feng's symbols "
                          f"({absent[:5]})")
        set_idx[name] = [uidx[g] for g in found]
        ctrl_idx[name] = np.array(sorted(uidx[g] for g in controls), dtype=np.int64)
        shown[name] = {"found": found, "missing": [g for g in s["genes"] if g not in uidx], "controls": controls}
    return set_idx, ctrl_idx, shown


def relabel_fixed(C: sparse.spmatrix, genes: Sequence[str], stage: np.ndarray, cfg: Mapping,
                  record_sets: Mapping) -> dict:
    """The E18 labelling on counts and stage labels only: each cell's set scores against row N54's control sets
    (`feng_relabel.set_score`), the yaml's assignment (`feng_relabel.assign`, its min_score, min_margin and split), and
    the yaml's checks on the E18 stage (the per-stage unassigned cap, the groups' cells, the fallback rule reported)."""
    stage = np.asarray(stage, dtype=str)
    if not np.isin(stage, list(STAGES)).all():
        raise Refused(f"a cell outside {list(STAGES)} was handed to the E18 labeller")
    L = fr.log1p_cp10k(C)
    set_idx, ctrl, shown = fixed_sets(genes, cfg, record_sets)
    scores = {name: fr.set_score(L, set_idx[name], ctrl[name]) for name in cfg["sets"]}
    labels = fr.assign(scores, cfg)
    stages = list(STAGES)
    low = np.vstack([scores[n] for n in fr.classes(cfg)]).max(0) < float(cfg["assign"]["min_score"])
    shares = fr.unassigned_shares(labels, stage, stages, cfg["assign"]["unassigned"])
    cap = fr.cap_check(shares, float(cfg["unassigned_cap"]))
    cells = fr.group_cells(labels, stage, stages, cfg["groups"])
    unassigned = cfg["assign"]["unassigned"]
    kinds = {s: {"no_class_reached_min_score": int(np.sum((stage == s) & low)),
                 "margin_failed": int(np.sum((stage == s) & (labels == unassigned) & ~low))} for s in stages}
    return {"labels": labels, "scores": scores, "unassigned_share": shares, "cap": cap, "group_cells": cells,
            "unassigned_kinds": kinds, "sets": shown,
            "label_counts": {s: {lab: int(np.sum((stage == s) & (labels == lab))) for lab in fr.all_labels(cfg)}
                             for s in stages}}


# ------------------------------------------------------------------------------------------------ the copied rows


def load_raw(path: Path) -> dict:
    """A statistics file's arrays as stored (no conversion), every field of STAT_FIELDS required."""
    with np.load(path, allow_pickle=False) as z:
        missing = [f for f in STAT_FIELDS if f not in z.files]
        if missing:
            raise Refused(f"{Path(path).name} lacks the fields {missing}")
        return {f: z[f] for f in STAT_FIELDS}


def combine(n54: Mapping, e18: Mapping) -> dict:
    """Row N54's rows first, unchanged, then the E18 rows, in one statistics table. Refused unless both share one gene
    list, N54's keys are all of COPIED_STAGES and the E18 keys all of STAGE."""
    if len(e18["keys"]) == 0:
        raise Refused(f"the E18 statistics hold no row: no {STAGE} cell was labelled")
    if [str(g) for g in n54["genes"]] != [str(g) for g in e18["genes"]]:
        raise Refused("the E18 statistics and row N54's are not on one gene list")
    if {str(k).split("|")[0] for k in n54["keys"]} - set(COPIED_STAGES):
        raise Refused(f"row N54's statistics hold a stage outside {list(COPIED_STAGES)}")
    if {str(k).split("|")[0] for k in e18["keys"]} - set(STAGES):
        raise Refused(f"the E18 statistics hold a stage other than {STAGE}")
    out = {"keys": np.array([*map(str, n54["keys"]), *map(str, e18["keys"])]), "genes": np.asarray(n54["genes"])}
    for f in ROW_FIELDS:
        a, b = np.asarray(n54[f]), np.asarray(e18[f])
        if a.dtype != b.dtype:
            b = b.astype(a.dtype)
        out[f] = np.concatenate([a, b], axis=0)
    return out


def copied_rows_identical(combined: Mapping, n54: Mapping) -> dict:
    """The combined table's E9.5 and E14.5 rows against row N54's file: the same gene list; exactly N54's keys at those
    stages, each once; every row field byte for byte (same dtype); and every other key at E18. Returns the checks and
    'identical' (all of them)."""
    keys = [str(k) for k in combined["keys"]]
    want = [str(k) for k in n54["keys"]]
    pos = {}
    for i, k in enumerate(keys):
        pos.setdefault(k, []).append(i)
    early = [k for k in keys if k.split("|")[0] in COPIED_STAGES]
    checks = {"genes": [str(g) for g in combined["genes"]] == [str(g) for g in n54["genes"]],
              "keys": sorted(early) == sorted(want) and all(len(pos.get(k, [])) == 1 for k in want),
              "other_stages": all(k.split("|")[0] in (*COPIED_STAGES, STAGE) for k in keys)}
    rows_ok = checks["keys"]
    if rows_ok:
        for f in ROW_FIELDS:
            a, b = np.asarray(combined[f]), np.asarray(n54[f])
            if a.dtype != b.dtype:
                rows_ok = False
                break
            for i, k in enumerate(want):
                if np.ascontiguousarray(a[pos[k][0]]).tobytes() != np.ascontiguousarray(b[i]).tobytes():
                    rows_ok = False
                    break
    checks["rows"] = bool(rows_ok)
    return {**checks, "identical": all(checks.values())}


# ------------------------------------------------------------------------------------------------ run


def run(cfg: Mapping, src: Path, dst: Path, gate: Mapping, n54_raw: bytes | None, n54_path: Path) -> dict:
    """Read N54's statistics (sha256 checked) and the clean reduction's E18 cells (manifest-checked, the files N54
    read), label them with N54's control sets, write the combined statistics (N54's rows copied, then E18's), the E18
    cells table and the manifest. Returns the manifest."""
    if not gate or gate.get("clean") is not True or not N_ROW.fullmatch(str(gate.get("row"))):
        raise Refused("no clean gate of an N row: the E18 labelling runs only through main(), after preflight")
    rec = n54_record(n54_raw, cfg)
    if fio.sha256_file(n54_path) != N54_STATS_SHA256:
        raise Refused(f"row N54's statistics file's sha256 is not {N54_STATS_SHA256}")
    n54 = load_raw(n54_path)
    inp = cfg["inputs"]
    paths = {k: src / inp[k] for k in ("counts", "genes", "obs")}
    entries = {k: fio.verify_manifest_output(p) for k, p in paths.items()}
    read_by_n54 = {k: ((rec.get("inputs") or {}).get(k) or {}).get("sha256") for k in paths}
    if any(entries[k]["sha256"] != read_by_n54[k] for k in paths):
        raise Refused("the clean reduction's files are not the ones row N54 read (their sha256 differ from its record)")
    obs = fr.read_obs(paths["obs"], cfg["obs_columns"])
    if not np.all(obs["strain"] == cfg["expect_strain"]):
        raise Refused(f"a cell's strain is not {cfg['expect_strain']}")
    genes = paths["genes"].read_text().splitlines()
    C = sparse.load_npz(paths["counts"]).tocsr()
    if C.shape != (len(obs["cell"]), len(genes)):
        raise Refused(f"counts {C.shape} do not match the obs table and genes")
    keep = np.flatnonzero(np.isin(obs["stage"], list(STAGES)))
    if keep.size == 0:
        raise Refused(f"the clean reduction holds no cell at {STAGE}")
    names, stage, zone = obs["cell"][keep], obs["stage"][keep], obs["zone"][keep]
    C = C[keep]
    res = relabel_fixed(C, genes, stage, cfg, rec["sets"])
    e18 = feng2022_geo.group_stats(C, stage, zone, res["labels"], names, genes)
    combined = combine(n54, e18)
    same = copied_rows_identical(combined, n54)
    if not same["identical"]:
        raise Refused(f"the copied E9.5 and E14.5 rows are not row N54's byte for byte: {same}")
    dst.mkdir(parents=True, exist_ok=True)
    fio.save_npz_atomic(dst / STATS, **combined)
    fr._write_cells(dst / CELLS, names, stage, zone, res["labels"], res["scores"])
    again = copied_rows_identical(load_raw(dst / STATS), n54)  # the file as written, read back
    if not again["identical"] or fio.sha256_file(n54_path) != N54_STATS_SHA256:
        raise Refused(f"the written file's copied rows or row N54's file changed: {again}")
    usable = bool(res["cap"]["passed"])
    outputs = {STATS: {"sha256": fio.sha256_file(dst / STATS), "usable": usable, "groups": int(len(combined["keys"])),
                       "e18_rows": int(len(e18["keys"])), "half_seed": feng2022_geo.HALF_SEED},
               CELLS: {"sha256": fio.sha256_file(dst / CELLS)}}
    manifest = {"source": SOURCE, "method": cfg["method"], "register_row": "X002",
                "config": {"file": fr._rel(fr.CONFIG), "sha256": fio.sha256_file(fr.CONFIG),
                           "status": cfg.get("status")},
                "stages": list(STAGES), "stages_copied": list(COPIED_STAGES),
                "controls_from": {"record": fr._rel(fr.RECORD), "sha256": gate.get("n54_record_sha256"),
                                  "rule": "row N54's recorded control genes per set, unchanged; none drawn here"},
                "n54_stats": {"file": fr._rel(n54_path), "sha256": N54_STATS_SHA256},
                "copied_rows": again,
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


def write_stop(dst: Path, cfg: Mapping, gate: Mapping, reason: str) -> dict:
    """A run that stopped after preflight: a manifest listing no output (no reader can take a file from it), usable
    false, and its copy at RECORD."""
    manifest = {"source": SOURCE, "method": cfg.get("method"), "register_row": "X002",
                "config": {"file": fr._rel(fr.CONFIG), "sha256": fio.sha256_file(fr.CONFIG)},
                "stages": list(STAGES), "stages_copied": list(COPIED_STAGES), "gate": dict(gate), "stopped": True,
                "reason": str(reason)[:500], "usable": False, "authors_labels_read": False, "outputs": {},
                "provenance": fio.provenance(CODE)}
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
        dst = fr.output_dir_ok(fio.REPO / E18_DIR)
        n54_path = n54_stats_path()
        gate = preflight(a.row, dst)
    except (Refused, fr.Refused) as e:
        raise SystemExit(f"refused: {e}") from None
    dst.mkdir(parents=True)
    (dst / fr.STARTED).write_text(json.dumps({"gate": gate}, indent=1) + "\n")
    try:
        man = run(cfg, src, dst, gate, fr.committed_at_head(fr.RECORD), n54_path)
    except BaseException as e:  # any stop after preflight is recorded; nothing is rerun under the row
        write_stop(dst, cfg, gate, f"{type(e).__name__}: {e}")
        print(f"E18 labelling stopped after preflight; commit {RECORD.relative_to(fio.REPO)} now")
        if not isinstance(e, Exception):
            raise
        return STOP_EXIT
    rec = write_record(dst / "MANIFEST.json")
    print(f"E18 labelling written; usable {man['usable']}; commit {rec.relative_to(fio.REPO)} now")
    return 0


if __name__ == "__main__":
    sys.exit(main())
