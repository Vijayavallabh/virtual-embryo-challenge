"""File WBH, the heart's world-blur hedge (row R69 of experiments/harness/TRIES.md, the captain's order of 2026-11-25,
ruling, no run): the banked half re-measured in 32 antithetic, spatially blocked worlds (vec.world_blur, with every fix
of HEART_MOON2_challenge.md).

  world_blur.py IN OUT --prior same_stage --blocks 32 --seed 1000001 --captain-ruling R69 [--data DIR]

The constants are fixed: no option changes them, and any other value is refused. Refused before anything is read
beyond the ledger and IN's bytes: a ruling other than R69; a ledger without exactly one row R69 whose first cell is the
captain's order ("ruling, no run") naming File WBH, its gate row H37 and `--captain-ruling R69`; a row H37 whose
outcome cell does not hold the gate's PASS line (experiments/harness/h37_row.py; the file is built only after its gate
row reads PASS); an IN not named T2heart_organ_ts125_tk_half.h5ad or not of the sha256 that
experiments/builds/results/hearth_8a8bc3a_h2.SHA256SUMS lists under that name (OTK125's half); OUT equal to IN, or an
OUT that exists. Then the prior from row R69's five files under DIR (default data/raw), the blur, and a refusal if obs
names, variable order or the bytes of obsm["spatial_3D"] differ from IN's. The file is written as halve.py writes
(anndata's write_h5ad, gzip) through vec.io.write_submission, which stamps it. The tool prints "ok" and the output's
sha256 only, never a cell count, and writes no record beside the file.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np

from vec import harness
from vec import world_blur as wb
from vec.boards import load_boards
from vec.io import write_submission

ROOT = Path(__file__).resolve().parents[2]
TRIES_FILE = ROOT / "experiments" / "harness" / "TRIES.md"
SUMS = ROOT / "experiments" / "builds" / "results" / "hearth_8a8bc3a_h2.SHA256SUMS"
IN_NAME = "T2heart_organ_ts125_tk_half.h5ad"
RULING = "R69"
PRIOR = "same_stage"
BLOCKS = wb.N_BLOCKS
SEED = harness.BUILD_SEED_OFFSET + 1
BOARD_KEY = "T2:heart:test_extrap"
COORD = "spatial_3D"


def _rows():
    spec = importlib.util.spec_from_file_location("h37_row", ROOT / "experiments" / "harness" / "h37_row.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256_of(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def recorded_sha256(sums: Path = SUMS, name: str = IN_NAME) -> str | None:
    """The one sha256 the record lists under `name` (every line naming it must agree); None otherwise."""
    try:
        lines = sums.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    found = {ln.split()[0] for ln in lines if len(ln.split()) == 2 and Path(ln.split()[1]).name == name}
    return found.pop() if len(found) == 1 else None


def refusals(a, tries: Path | None = None, sums: Path | None = None) -> None:
    """Every refusal made before the prior or the blur: SystemExit with the reason."""
    if a.captain_ruling != RULING:
        sys.exit(f"world_blur.py is refused: it runs only with --captain-ruling {RULING} (file WBH's ruling row)")
    if a.prior != PRIOR or a.blocks != BLOCKS or a.seed != SEED:
        sys.exit(f"world_blur.py is refused: the constants are fixed (--prior {PRIOR} --blocks {BLOCKS} --seed {SEED})")
    rows = _rows()
    try:
        text = Path(tries or TRIES_FILE).read_text(encoding="utf-8")
    except OSError:
        text = ""
    if not rows.ruling_ok(text):
        sys.exit(f"world_blur.py is refused: TRIES.md does not hold exactly one row {RULING} whose first cell is the "
                 "captain's order (\"ruling, no run\") naming File WBH, its gate row H37 and "
                 f"`--captain-ruling {RULING}`")
    if (rows.outcome(text) or ("OTHER", None))[0] != "PASS":
        sys.exit("world_blur.py is refused: row H37's outcome cell does not hold the gate's PASS line (the file is "
                 "built only after its gate row reads PASS)")
    src, out = Path(a.src), Path(a.out)
    want = recorded_sha256(Path(sums or SUMS))
    if src.name != IN_NAME or not src.is_file() or want is None or sha256_of(src) != want:
        sys.exit(f"world_blur.py is refused: IN must be {IN_NAME} of the sha256 {SUMS.name} lists under that name")
    if os.path.realpath(src) == os.path.realpath(out):
        sys.exit("world_blur.py is refused: OUT must not be IN")
    if out.exists():
        sys.exit(f"world_blur.py is refused: {out} exists")


def blurred_file(src: Path, s: np.ndarray, genes: list[str]) -> ad.AnnData:
    """IN with its X re-measured at the file's seed; refused if IN's variables are not the board's genes in its order
    (the prior's order) or if anything but X would change."""
    a = ad.read_h5ad(src)
    if list(a.var_names) != list(genes):
        sys.exit("world_blur.py is refused: IN's variables are not the board's genes in its order")
    b = a.copy()
    b.X = wb.blur(a.X, np.asarray(a.obsm[COORD])[:, :3], s, SEED).X
    same = (list(b.obs_names) == list(a.obs_names) and list(b.var_names) == list(a.var_names)
            and np.asarray(b.obsm[COORD]).tobytes() == np.asarray(a.obsm[COORD]).tobytes()
            and list(b.obs.columns) == list(a.obs.columns) and b.X.dtype == a.X.dtype and type(b.X) is type(a.X))
    if not same:
        sys.exit("world_blur.py is refused: obs names, variable order, coordinates or X's type would change")
    return b


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src", metavar="IN")
    ap.add_argument("out", metavar="OUT")
    ap.add_argument("--prior", required=True)
    ap.add_argument("--blocks", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--captain-ruling", required=True)
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    return ap


def main(argv=None, tries: Path | None = None, sums: Path | None = None, genes=None) -> str:
    a = build_parser().parse_args(argv)
    refusals(a, tries, sums)
    genes = list(genes) if genes is not None else list(load_boards()[BOARD_KEY].genes)
    prior = wb.load_prior(Path(a.data), genes)
    if not np.isfinite(prior.s).all():
        sys.exit("world_blur.py is refused: the prior cannot be formed from row R69's files")
    b = blurred_file(Path(a.src), prior.s, genes)
    write_submission(b, a.out)
    digest = sha256_of(a.out)
    print("ok")
    print(digest)
    return digest


if __name__ == "__main__":
    main(sys.argv[1:])
