"""JAX mouse time-lapse (Qiu et al. 2024): window filter by embryo ID.

Keep rule (docs/DATA_POLICY.md; docs/PLAN.md section 4.1): embryos 1-22 (labelled E8.5-E9.5, 0-24 somites) and 54-74
(E14.25-P0). Everything else is dropped before anything is written:
  23-24   boundary: staged E9.5 by somites (25-26) but harvested in a litter labelled E9.75
  25-49   inside the T1 / T2-heart window (litters E9.75 to E13.5); never downloaded
  50-53   boundary: labelled E13.75-E14.0, just after the window
  75+     not part of the legal late anchor (includes the JAX validation embryos)
Embryos 1-12 are T1-only: their E8.5 bin (0-12 somites) overlaps the T2-heart interpolation window (E8.25, E8.75), and
embryo 1 (0 somites) may sit in the T2-embryo window (E7.25, E8.0). Naming any task besides T1 with --tasks drops them.
Every kept cell row is also re-checked: its day label must lie outside the windows of every named task (abort if not).

Inputs: the legal subsets already on the DGX (paths in sources.yaml), written by the 2026-10-18 research run:
  cell tables   df_cell.legalT1.csv (embryos 1-22, 50-74), lpm_legalT1.csv (LPM sub-labels, same embryos),
                jax_labels.csv (refined labels, embryos 1-22, 54-74): streamed row by row on the embryo_id column
  gene table    df_gene.csv (gene IDs, symbols, chromosome; no cells): copied
  per embryo    out/jax/embryo_N.npz (per-label sufficient statistics) and ext/jax/embryo/embryo_N_gene_count.rds (raw
                counts): gated by the embryo number in the file name (and, for the npz, the embryo field inside it);
                files that pass are hard-linked into the destination (copied if on another filesystem)

Usage (DGX, repo root; a few minutes, mostly hashing 25 GB of count matrices):
  uv run python experiments/external/jax_qiu2024.py [--src DIR] [--dst data/external/jax_qiu2024] [--tasks T1]
      [--skip-rds]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import filterio as fio  # noqa: E402
import windows  # noqa: E402

SOURCE = "jax_qiu2024"
KEEP_EMBRYOS = frozenset(range(1, 23)) | frozenset(range(54, 75))
T1_ONLY_EMBRYOS = frozenset(range(1, 13))
BOUNDARY_EMBRYOS = frozenset({23, 24, 50, 51, 52, 53})
IN_WINDOW_EMBRYOS = frozenset(range(25, 50))
CODE = [Path(__file__), Path(windows.__file__), Path(fio.__file__)]
_EMBRYO = re.compile(r"^embryo_(\d+)$")


def embryo_number(embryo_id: str) -> int:
    m = _EMBRYO.match(str(embryo_id).strip())
    if not m:
        raise ValueError(f"malformed embryo_id {embryo_id!r}")
    return int(m.group(1))


def embryos_for(tasks=("T1",)) -> frozenset[int]:
    """Embryos usable for every task in `tasks`: the keep set, minus the T1-only embryos if any other task is named."""
    tasks = windows.check_tasks(tasks)
    return KEEP_EMBRYOS if set(tasks) == {"T1"} else KEEP_EMBRYOS - T1_ONLY_EMBRYOS


def keep_embryo(n: int, tasks=("T1",)) -> tuple[bool, str]:
    if n in embryos_for(tasks):
        return True, "kept"
    if n in T1_ONLY_EMBRYOS:
        return False, "dropped: T1-only embryo (1-12) for a non-T1 task"
    if n in IN_WINDOW_EMBRYOS:
        return False, "dropped: in-window embryo (25-49)"
    if n in BOUNDARY_EMBRYOS:
        return False, "dropped: boundary embryo (23-24, 50-53)"
    return False, "dropped: outside the keep set"


def decide(values: tuple[str, ...], tasks=("T1",)) -> tuple[bool, str]:
    """Keep rule for one cell row, from its (embryo_id, day) values."""
    embryo_id, day = values
    keep, reason = keep_embryo(embryo_number(embryo_id), tasks)
    if keep and not windows.allowed_label(day, tasks):  # the embryo list already excludes every window
        raise AssertionError(f"{embryo_id} is in the keep set but its day label {day!r} lies in a window of {tasks}")
    return keep, reason


def filter_cells(src: str | Path, dst: str | Path, tasks=("T1",)) -> dict:
    return fio.filter_csv(src, dst, ["embryo_id", "day"], lambda v: decide(v, tasks))


def gate_files(src_pattern: Path, dst_pattern: Path, tasks=("T1",), check_npz: bool = False,
               overwrite: bool = False, hash_files: bool = True) -> dict:
    """Place every per-embryo file whose embryo passes the keep rule; files of other embryos are never touched.

    `src_pattern` and `dst_pattern` contain '{n}'. Only embryos that pass are looked up, so a file of a dropped embryo
    in the source directory is not opened, copied or hashed (it is only listed, by name, under `skipped`).
    """
    keep = sorted(embryos_for(tasks))
    rx = re.compile("^" + re.escape(src_pattern.name).replace(r"\{n\}", r"(\d+)") + "$")
    present = {int(m.group(1)): p for p in src_pattern.parent.iterdir() if (m := rx.match(p.name))}
    skipped = sorted(n for n in present if n not in keep)
    missing = [n for n in keep if n not in present]
    placed = {}
    for n in keep:
        if n not in present:
            continue
        src = present[n]
        if check_npz:
            z = np.load(src, allow_pickle=True)  # written by the research run; check the embryo it holds
            inner = str(z["embryo"])
            if inner != f"embryo_{n}":
                raise ValueError(f"{src}: file name says embryo {n}, contents say {inner}")
        dst = Path(str(dst_pattern).format(n=n))
        mode = fio.link_or_copy(src, dst, overwrite=overwrite)
        placed[n] = {"src": str(src), "dst": str(dst), "mode": mode,
                     "sha256": fio.sha256_file(dst) if hash_files else None, "bytes": dst.stat().st_size}
    return {"kept": sorted(placed), "missing_from_source": missing, "skipped_not_kept": skipped, "files": placed}


def run(src: Path, dst: Path, cfg: dict, tasks=("T1",), skip_rds: bool = False, overwrite: bool = False) -> dict:
    dst.mkdir(parents=True, exist_ok=True)
    inputs, outputs = {}, {}
    for rel_src, rel_dst in cfg["cell_tables"].items():
        p = src / rel_src
        rep = filter_cells(p, dst / rel_dst, tasks)
        inputs[rel_src] = {"path": str(p), "sha256": fio.sha256_file(p)}
        outputs[rel_dst] = {"from": rel_src, "sha256": fio.sha256_file(dst / rel_dst), **rep}
    for rel_src, rel_dst in cfg["gene_tables"].items():
        p = src / rel_src
        mode = fio.link_or_copy(p, dst / rel_dst, overwrite=overwrite)
        inputs[rel_src] = {"path": str(p), "sha256": fio.sha256_file(p)}
        outputs[rel_dst] = {"from": rel_src, "mode": mode, "sha256": inputs[rel_src]["sha256"], "note": "no cell rows"}
    per_embryo = {}
    for rel_src, rel_dst in cfg["per_embryo"].items():
        if skip_rds and rel_src.endswith(".rds"):
            per_embryo[rel_src] = "skipped (--skip-rds)"
            continue
        rep = gate_files(src / rel_src, dst / rel_dst, tasks, check_npz=rel_src.endswith(".npz"), overwrite=overwrite)
        per_embryo[rel_src] = {k: v for k, v in rep.items() if k != "files"}
        for f in rep["files"].values():
            outputs[str(Path(f["dst"]).relative_to(dst))] = {"from": f["src"], "mode": f["mode"], "sha256": f["sha256"],
                                                             "bytes": f["bytes"]}
    return {"inputs": inputs, "per_embryo": per_embryo, "outputs": outputs}


def main(argv=None) -> int:
    cfg = yaml.safe_load((fio.HERE / "sources.yaml").read_text())[SOURCE]
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src", default=cfg["src_dir"], help="root of the legal subsets")
    ap.add_argument("--dst", default=cfg["dst"], help="output directory (must be under data/, which is gitignored)")
    ap.add_argument("--tasks", nargs="+", default=list(cfg["tasks"]), help="tasks the data may serve (union)")
    ap.add_argument("--skip-rds", action="store_true", help="leave out the raw per-embryo count matrices (25 GB)")
    ap.add_argument("--overwrite", action="store_true", help="replace destination files that differ from the source")
    a = ap.parse_args(argv)
    tasks = windows.check_tasks(a.tasks)
    src = Path(a.src).expanduser()
    dst = (fio.REPO / a.dst).resolve() if not Path(a.dst).is_absolute() else Path(a.dst)
    if fio.REPO / "data" not in dst.parents:
        raise SystemExit(f"--dst must be inside {fio.REPO / 'data'} (gitignored); got {dst}")
    res = run(src, dst, cfg["inputs"], tasks, a.skip_rds, a.overwrite)
    manifest = {"source": SOURCE, "citation": cfg["citation"], "url": cfg["url"], "tasks": list(tasks),
                "keep_rule": cfg["keep"], "keep_embryos": sorted(embryos_for(tasks)),
                "t1_only_embryos": sorted(T1_ONLY_EMBRYOS), **res, "provenance": fio.provenance(CODE)}
    path = fio.write_manifest(dst, manifest)
    for name, o in res["outputs"].items():
        if "rows_in" in o:
            print(name, {k: o[k] for k in ("rows_in", "rows_kept", "rows_dropped", "by_reason")})
    for name, o in res["per_embryo"].items():
        print(name, o if isinstance(o, str) else {k: (len(v) if k == "kept" else v) for k, v in o.items()})
    print("manifest:", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
