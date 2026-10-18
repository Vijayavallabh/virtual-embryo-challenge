"""Preflight for run_p3.sh (DGX): stop before any measurement unless every input is usable.

Checks, in order (exit code in brackets):
  [2] every input file exists;
  [3] every file has obs["celltype"] and every label in it is in configs/lineage_map.yaml (the mapping is
      vec.prep.add_lineage itself, applied to obs only, so it is fast); new labels must be mapped by the rule in the
      map's header and committed first;
  [4] the files that define the procedure have no uncommitted changes (pre-registration): configs/lineage_map.yaml,
      experiments/coords/size_procedure.py, src/vec/prep.py, src/vec/spatial.py (skip with --allow-uncommitted; the
      size JSON records uncommitted changes either way).
Prints one line per file (labels, core cells under the map). Usage:
  uv run python experiments/coords/p3_preflight.py [--allow-uncommitted] FILE...
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PROCEDURE_FILES = ("configs/lineage_map.yaml", "experiments/coords/size_procedure.py", "src/vec/prep.py",
                   "src/vec/spatial.py", "experiments/coords/handedness.py")


def label_check(path: Path) -> tuple[list[str], str]:
    """(unmapped labels, one-line summary); raises KeyError when the file has no celltype column."""
    import anndata as ad

    from vec.prep import add_lineage

    a = ad.read_h5ad(path, backed="r")
    try:
        if "celltype" not in a.obs:
            raise KeyError(f"{path}: no obs['celltype'] column")
        b = ad.AnnData(obs=a.obs[["celltype"]].copy())
    finally:
        a.file.close()
    unmapped = add_lineage(b, source=path.name.removesuffix(".gz").removesuffix(".h5ad"))
    core = b.obs["core"].astype(bool)
    core_labels = b.obs.loc[core, "celltype"].astype(str).value_counts()
    return unmapped, (f"{path}: {b.n_obs} cells, {b.obs['celltype'].nunique()} labels, core {int(core.sum())} cells "
                      f"({', '.join(f'{k} {v}' for k, v in core_labels.items())})")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("--allow-uncommitted", action="store_true")
    a = ap.parse_args(argv)
    missing = [str(f) for f in a.files if not f.is_file()]
    if missing:
        print("PREFLIGHT FAILED: missing input files:\n  " + "\n  ".join(missing), file=sys.stderr)
        return 2
    bad = {}
    for f in a.files:
        try:
            unmapped, line = label_check(f)
        except KeyError as e:
            print(f"PREFLIGHT FAILED: {e}", file=sys.stderr)
            return 3
        print(line)
        if unmapped:
            bad[str(f)] = unmapped
    if bad:
        print("PREFLIGHT FAILED: labels missing from configs/lineage_map.yaml (map them by the rule in its header, "
              "commit, sync, re-run):", file=sys.stderr)
        for f, labs in bad.items():
            print(f"  {f}: {labs}", file=sys.stderr)
        return 3
    dirty = subprocess.run(["git", "status", "--porcelain", "--", *PROCEDURE_FILES], cwd=REPO, capture_output=True,
                           text=True, check=False).stdout.rstrip("\n")
    if dirty:
        msg = f"uncommitted changes to the files that define the procedure:\n{dirty}"
        if not a.allow_uncommitted:
            print(f"PREFLIGHT FAILED: {msg}\n  commit them first (pre-registration), or pass --allow-uncommitted",
                  file=sys.stderr)
            return 4
        print(f"WARNING (allowed): {msg}")
    print("preflight ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
