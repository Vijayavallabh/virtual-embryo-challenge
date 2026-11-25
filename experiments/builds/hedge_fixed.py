#!/usr/bin/env python
"""Row R68's four hedge files (the captain's order of 2026-11-25: "Yes to the four hedge files of HEDGE_FIXED_design.md
as named"), as that score-blind design's section 3 fixes them: each row of a file comes whole from one of two parents
fixed by recipe, 50:50, on the rows, names and positions of an anchor file.

  file                      board                  anchor  first  second  pick
  T3_c2_s2_mxe2.h5ad        T3:bcatenin            C2      C2     xe2     parity
  T3_s2_mx_xe2_ck.h5ad      T3:bcatenin            C2      xe2    ck      parity
  T2embryo_s2x_ovmxso.h5ad  T2:embryo:test_interp  O       O      SO      parity
  T2embryo_s2x_ovbkso.h5ad  T2:embryo:test_interp  O       O      SO      blocks

The parents by the sha256 their records' committed lists hold (experiments/builds/results/<record>.SHA256SUMS): C2
c2-s2.h5ad (t3two_fb01487), xe2 xall-s2.h5ad (t3k23_49f3fd5), ck ck-s2.h5ad (t3ck_20a0e1a), O O/T2embryo_s2x_ov.h5ad
and SO SO/T2embryo_s2xs_ov.h5ad (embovsb_20b2381_w1). The picks; no count enters either:

  parity  row i from the first parent where i is even, from the second where i is odd
  blocks  the anchor's (O's) coordinates centred, put on their principal axes (each axis's sign set so that its largest
          loading is positive) and divided by their RMS radius (vec.t2_build.rms_radius): row i from the second parent
          where floor(z1) + floor(z2) + floor(z3) is odd, else from the first

Row i of X is the picked parent's row i, bit for bit; obs and obsm are the anchor's (names cell_0 ... cell_<n-1>,
obsm["spatial_3D"] float32), so count, order and positions are the anchor's. The board's own writer writes the file:
T3 vec.io make_submission, write_submission (stamped) and check_submission; the embryo vec.t2_build.build_submission
(stamped and checked). The written file is read back and must hold the anchor's obs and obsm and the picked rows. The
sidecar (the file's name with .json for .h5ad) holds the names, the sha256s, the pick, this script's sha256 and the
ruling row, never a count; the script prints "ok" and the file's sha256, nothing else.

Exit status 2, nothing written (refused): no --captain-ruling R68 (or another row); a ledger
(experiments/harness/TRIES.md) without exactly one row R68 whose first cell passes
build_outside_form.captain_order_names for the file's name in backticks and holds the name, this script's path and
--captain-ruling R68 as whole backticked strings; a name not in the table; an --out not named as the file, or one
that exists (or its sidecar or stamp); a tree that is not clean (vec.provenance's rule); an anchor or parent that is
not a file with its pinned sha256, or the parents in the wrong order; a parent whose n_obs, obs_names, var_names or
obsm["spatial_3D"] is not the anchor's (np.array_equal, dtype included); var_names other than the board's genes; a
parent's X that is not dense float32, finite and >= 0; an anchor whose obs or obsm the board's writer would not give
back (other obs names, obs columns, another obsm key, coordinates that are not (n, 3) float32, a cloud without
spread). Exit status 1: the writer's own refusal, a failed format check or a file that does not read back as built
(messages with digits masked).

  uv run --frozen --no-sync python -I -B experiments/builds/hedge_fixed.py build --captain-ruling R68 \\
      --name T3_c2_s2_mxe2.h5ad --anchor <t3two>/c2-s2.h5ad --first <t3two>/c2-s2.h5ad \\
      --second <t3k23>/xall-s2.h5ad --out runs/hedgefix_<c7>_<RUN>/a/T3_c2_s2_mxe2.h5ad

experiments/builds/hedge_fixed_line.sh builds each file twice on the DGX and compares the two.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # this checkout
sys.path.insert(0, str(REPO / "src"))  # python -I: this checkout's vec, never the one a shared environment names
import vec.provenance  # noqa: E402, F401, I001  (first: records this build's start before any method code)

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
from typing import NamedTuple  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402

from vec import provenance  # noqa: E402
from vec.boards import BOARDS_DIR, Board, get_board  # noqa: E402
from vec.io import COORD_KEY, check_submission, make_submission, write_submission  # noqa: E402
from vec.t2_build import build_submission, rms_radius  # noqa: E402

SCRIPT = "experiments/builds/hedge_fixed.py"
SCRIPT_PATH = Path(__file__).resolve()
TRIES_FILE = REPO / "experiments" / "harness" / "TRIES.md"
OUTSIDE_FORM_DIR = REPO / "experiments" / "t3_decomp"  # build_outside_form.py, whose captain_order_names is reused
RULING_ROW = "R68"
DESIGN = "HEDGE_FIXED_design.md (2026-11-25), section 3"
EXIT_FAILED, EXIT_REFUSED = 1, 2
T3_BOARD, EMBRYO_BOARD = "T3:bcatenin", "T2:embryo:test_interp"
C2_SHA = "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"  # t3two_fb01487.SHA256SUMS
XE2_SHA = "77ff0b947f7f409c9ce171bdc59a44de3d3c2ef2f3e5a7332011a459ce9ee9c3"  # t3k23_49f3fd5.SHA256SUMS
CK_SHA = "e470aaeeab28918e157fb57cdc771122af35e46f0916cd1c6d63650c31faae26"  # t3ck_20a0e1a.SHA256SUMS
O_SHA = "1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea"  # embovsb_20b2381_w1.SHA256SUMS
SO_SHA = "77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165"  # embovsb_20b2381_w1.SHA256SUMS
PARENTS = {  # sha256 -> (the design's name, the file as its record lists it, the record)
    C2_SHA: ("C2", "c2-s2.h5ad", "t3two_fb01487"),
    XE2_SHA: ("xe2", "xall-s2.h5ad", "t3k23_49f3fd5"),
    CK_SHA: ("ck", "ck-s2.h5ad", "t3ck_20a0e1a"),
    O_SHA: ("O", "O/T2embryo_s2x_ov.h5ad", "embovsb_20b2381_w1"),
    SO_SHA: ("SO", "SO/T2embryo_s2xs_ov.h5ad", "embovsb_20b2381_w1"),
}


class Hedge(NamedTuple):
    board: str
    parents: tuple[str, str]  # the first and the second parent's sha256, in that order
    pick: str  # "parity" or "blocks"
    anchor: str  # the anchor's sha256: its rows, names and positions are the file's


TABLE = {  # the pinned table: name -> board, parents in order, pick, anchor
    "T3_c2_s2_mxe2.h5ad": Hedge(T3_BOARD, (C2_SHA, XE2_SHA), "parity", C2_SHA),
    "T3_s2_mx_xe2_ck.h5ad": Hedge(T3_BOARD, (XE2_SHA, CK_SHA), "parity", C2_SHA),
    "T2embryo_s2x_ovmxso.h5ad": Hedge(EMBRYO_BOARD, (O_SHA, SO_SHA), "parity", O_SHA),
    "T2embryo_s2x_ovbkso.h5ad": Hedge(EMBRYO_BOARD, (O_SHA, SO_SHA), "blocks", O_SHA),
}
PICK_RULES = {
    "parity": "row i from the first parent where i is even, from the second where i is odd",
    "blocks": ("the anchor's coordinates centred, put on their principal axes (each axis's sign set so that its "
               "largest loading is positive) and divided by their RMS radius: row i from the second parent where "
               "floor(z1) + floor(z2) + floor(z3) is odd, else from the first"),
}
WRITERS = {T3_BOARD: "vec.io make_submission, write_submission, check_submission",
           EMBRYO_BOARD: "vec.t2_build.build_submission"}


class Refused(Exception):
    """Refused: exit status 2, nothing written."""


class Failed(Exception):
    """The writer's refusal, a failed format check or a file that does not read back as built: exit status 1."""


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _masked(msg: str) -> str:
    return re.sub(r"\d", "#", msg)


def hedge(name) -> Hedge:
    if name not in TABLE:
        raise Refused(f"--name {name!r} is not one of row {RULING_ROW}'s four files ({', '.join(TABLE)})")
    return TABLE[name]


def _outside_form():
    """build_outside_form, imported only here (its captain_order_names is reused, not copied)."""
    if str(OUTSIDE_FORM_DIR) not in sys.path:
        sys.path.append(str(OUTSIDE_FORM_DIR))
    import build_outside_form

    return build_outside_form


def first_cell(row: str, path: Path) -> str | None:
    """The first cell of the ledger's single row `row`, read as captain_order_names reads it; None for none or two."""
    if not path.is_file():
        return None
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    return lines[0][2:].split("|", 1)[0] if len(lines) == 1 else None


def names_whole(cell: str, s: str) -> bool:
    """Whether the cell holds `s` as one whole backtick-delimited string."""
    return s in re.findall(r"`([^`]*)`", cell)


def check_ruling(ruling, name: str, tries: Path | None = None) -> dict:
    """--captain-ruling R68 and the ledger's single row R68 naming the file; returns its record, else Refused."""
    if ruling is None:
        raise Refused(f"no --captain-ruling {RULING_ROW}: the files of this script are built only with it")
    if ruling != RULING_ROW:
        raise Refused(f"--captain-ruling {ruling!r} is refused: the files of this script are built only with "
                      f"--captain-ruling {RULING_ROW}")
    p = Path(TRIES_FILE if tries is None else tries)
    cell = first_cell(RULING_ROW, p)
    if not (_outside_form().captain_order_names(RULING_ROW, f"`{name}`", p) and cell is not None
            and all(names_whole(cell, s) for s in (name, SCRIPT, f"--captain-ruling {RULING_ROW}"))):
        raise Refused(f"experiments/harness/TRIES.md holds no single row {RULING_ROW} whose first cell is the "
                      f"captain's order (\"ruling, no run\") naming `{name}`, `{SCRIPT}` and "
                      f"`--captain-ruling {RULING_ROW}`")
    return {"row": RULING_ROW, "ledger": "experiments/harness/TRIES.md", "ledger_sha256": sha256(p)}


def check_out(out, name: str) -> Path:
    p = Path(out)
    if p.name != name:
        raise Refused(f"--out {out} is not named {name}")
    for q in (p, p.with_suffix(".json"), provenance.stamp_path(p)):
        if q.exists() or q.is_symlink():
            raise Refused(f"{q} exists: nothing is replaced")
    return p


def check_clean() -> dict:
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    if commit == "unknown" or unclean is None:
        raise Refused("git cannot tell the tree's state (the stamp would not say clean)")
    if unclean:
        raise Refused(f"the tree is not clean (e.g. {unclean[:3]}): commit first")
    return {"commit": commit, "clean": True}


def check_pins(h: Hedge, paths: dict[str, Path]) -> dict[str, str]:
    """Each of anchor, first and second a file with its pinned sha256, the parents in order; returns the sha256s."""
    for role, p in paths.items():
        if not p.is_file():
            raise Refused(f"--{role} {p} is not a file")
    got = {role: sha256(p) for role, p in paths.items()}
    if h.parents[0] != h.parents[1] and (got["first"], got["second"]) == h.parents[::-1]:
        raise Refused("the parents are in the wrong order: --first holds the second parent and --second the first")
    for role, want in (("anchor", h.anchor), ("first", h.parents[0]), ("second", h.parents[1])):
        if got[role] != want:
            raise Refused(f"--{role} {paths[role]}: sha256 {got[role][:12]} is not its pin {want[:12]} "
                          f"({PARENTS.get(want, ('?',))[0]})")
    return got


def _same(a, b) -> bool:
    """np.array_equal, dtype included."""
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and bool(np.array_equal(a, b))


def check_like_anchor(role: str, a: ad.AnnData, anchor: ad.AnnData) -> None:
    if a.n_obs != anchor.n_obs:
        raise Refused(f"--{role}: its n_obs is not the anchor's")
    if not _same(a.obs_names, anchor.obs_names):
        raise Refused(f"--{role}: its obs_names are not the anchor's (np.array_equal, dtype included)")
    if not _same(a.var_names, anchor.var_names):
        raise Refused(f"--{role}: its var_names are not the anchor's (np.array_equal, dtype included)")
    if COORD_KEY not in a.obsm or COORD_KEY not in anchor.obsm or not _same(a.obsm[COORD_KEY], anchor.obsm[COORD_KEY]):
        raise Refused(f"--{role}: its obsm['{COORD_KEY}'] is not the anchor's (np.array_equal, dtype included)")


def check_genes(role: str, a: ad.AnnData, board: Board) -> None:
    if [str(g) for g in a.var_names] != list(board.genes):
        raise Refused(f"--{role}: its var_names are not the genes of {board.key}, in their order")


def check_x(role: str, a: ad.AnnData) -> np.ndarray:
    X = a.X
    if sparse.issparse(X) or not isinstance(X, np.ndarray) or isinstance(X, np.matrix):
        raise Refused(f"--{role}: X is not a dense array")
    if X.dtype != np.float32:
        raise Refused(f"--{role}: X is {X.dtype}, not float32")
    if not np.isfinite(X).all():
        raise Refused(f"--{role}: X holds NaN or infinite values")
    if (X < 0).any():
        raise Refused(f"--{role}: X holds negative values")
    return X


def check_anchor_form(anchor: ad.AnnData) -> None:
    """The anchor's obs and obsm are what the board's writer gives back: names cell_0 ... cell_<n-1>, no obs column,
    obsm spatial_3D only, (n, 3) float32, with spread."""
    names = np.array([f"cell_{i}" for i in range(anchor.n_obs)], dtype=object)
    if not _same(anchor.obs_names, names):
        raise Refused("the anchor's obs names are not cell_0 ... cell_<n-1>, the board's writer's")
    if anchor.obs.shape[1]:
        raise Refused("the anchor's obs holds columns the board's writer does not write")
    if set(anchor.obsm.keys()) != {COORD_KEY}:
        raise Refused(f"the anchor's obsm holds other keys than {COORD_KEY}")
    C = anchor.obsm[COORD_KEY]
    if not isinstance(C, np.ndarray) or C.dtype != np.float32 or C.shape != (anchor.n_obs, 3):
        raise Refused(f"the anchor's obsm['{COORD_KEY}'] is not (n, 3) float32")
    if not (np.isfinite(C).all() and rms_radius(C) > 0):
        raise Refused(f"the anchor's obsm['{COORD_KEY}'] is not finite with spread")


def canonical(coords) -> np.ndarray:
    """Centred, on the principal axes (by decreasing variance; each axis's sign set so that its largest loading is
    positive), divided by the RMS radius. Sums and products in a fixed order, no threaded BLAS."""
    C = np.asarray(coords, dtype=np.float64)[:, :3]
    X = C - C.mean(0)
    cov = np.array([[float(np.sum(X[:, j] * X[:, k])) for k in range(3)] for j in range(3)])
    V = np.linalg.eigh(cov)[1][:, ::-1].copy()  # columns: the axes
    for k in range(3):
        if V[int(np.argmax(np.abs(V[:, k]))), k] < 0:
            V[:, k] = -V[:, k]
    Z = np.column_stack([X[:, 0] * V[0, k] + X[:, 1] * V[1, k] + X[:, 2] * V[2, k] for k in range(3)])
    return Z / rms_radius(C)


def second_rows(h: Hedge, anchor: ad.AnnData) -> np.ndarray:
    """True where row i comes from the second parent."""
    if h.pick == "parity":
        return np.arange(anchor.n_obs) % 2 == 1
    if h.pick == "blocks":
        return np.floor(canonical(anchor.obsm[COORD_KEY])).astype(np.int64).sum(1) % 2 == 1
    raise ValueError(f"unknown pick {h.pick!r}")


def write(h: Hedge, board: Board, X: np.ndarray, anchor: ad.AnnData, out: Path) -> None:
    """The board's own writer: the file and its stamp, checked; Failed when the writer refuses or the check fails."""
    coords = anchor.obsm[COORD_KEY]
    try:
        if h.board == T3_BOARD:
            write_submission(make_submission(X, board.genes, coords), out)
            rep = check_submission(out, board)
            if not rep.ok:
                raise Failed(f"the format check failed: {rep.errors}")
        elif h.board == EMBRYO_BOARD:
            build_submission(board, X, coords, genes=list(board.genes), out=out)
        else:
            raise Failed(f"no writer for {h.board}")
    except ValueError as e:
        raise Failed(f"the board's writer refused: {e}") from None


def check_written(out: Path, anchor: ad.AnnData, X: np.ndarray) -> None:
    b = ad.read_h5ad(out)
    if not (_same(b.obs_names, anchor.obs_names) and _same(b.var_names, anchor.var_names) and b.obs.shape[1] == 0
            and set(b.obsm.keys()) == {COORD_KEY} and _same(b.obsm[COORD_KEY], anchor.obsm[COORD_KEY])
            and isinstance(b.X, np.ndarray) and _same(b.X, X)):
        raise Failed("the written file does not hold the anchor's obs and obsm and the picked rows")


def _parent(path: Path, sha: str) -> dict:
    key, file, record = PARENTS[sha]
    return {"parent": key, "file": file, "record": f"experiments/builds/results/{record}.SHA256SUMS",
            "path": str(path), "sha256": sha}


def build(name, anchor, first, second, out, captain_ruling, tries: Path | None = None) -> dict:
    """One of the four files and its sidecar; returns the sidecar. Refused or Failed as the module says."""
    h = hedge(name)
    ruled = check_ruling(captain_ruling, name, tries)
    out_p = check_out(out, name)
    tree = check_clean()
    paths = {"anchor": Path(anchor), "first": Path(first), "second": Path(second)}
    got = check_pins(h, paths)
    board = get_board(h.board, BOARDS_DIR)
    read: dict[str, ad.AnnData] = {}  # by sha256: a file named twice is read once
    for role, p in paths.items():
        if got[role] not in read:
            read[got[role]] = ad.read_h5ad(p)
    files = {role: read[got[role]] for role in paths}
    A = files["anchor"]
    check_genes("anchor", A, board)
    check_anchor_form(A)
    X = {}
    for role in ("first", "second"):
        check_like_anchor(role, files[role], A)
        check_genes(role, files[role], board)
        X[role] = check_x(role, files[role])
    pick = second_rows(h, A)
    mixed = np.where(pick[:, None], X["second"], X["first"])
    write(h, board, mixed, A, out_p)
    check_written(out_p, A, mixed)
    side = {"script": SCRIPT, "script_sha256": sha256(SCRIPT_PATH), "design": DESIGN, "name": name,
            "board": h.board, "pick": h.pick, "pick_rule": PICK_RULES[h.pick],
            **{role: _parent(paths[role], got[role]) for role in ("anchor", "first", "second")},
            "obs_obsm": "the anchor's", "writer": WRITERS[h.board], "captain_ruling": ruled, "tree_at_start": tree,
            "output": {"file": name, "sha256": sha256(out_p), "format_check": "passed"}}
    out_p.with_suffix(".json").write_text(json.dumps(side, indent=1) + "\n")
    return side


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="hedge_fixed.py", description="Row R68's four hedge files (HEDGE_FIXED_design)")
    sub = ap.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="build one of the four files")
    b.add_argument("--captain-ruling", default=None, help=f"{RULING_ROW}, the only row that orders these files")
    b.add_argument("--name", required=True, help="one of the four files' names")
    b.add_argument("--anchor", required=True, help="the anchor (its rows, names and positions)")
    b.add_argument("--first", required=True, help="the first parent")
    b.add_argument("--second", required=True, help="the second parent")
    b.add_argument("--out", required=True, help="the file to write, named as --name (sidecar: .json)")
    a = ap.parse_args(argv)
    try:
        side = build(a.name, a.anchor, a.first, a.second, a.out, a.captain_ruling)
    except Refused as e:
        print(f"refused, nothing written: {e}", file=sys.stderr)
        return EXIT_REFUSED
    except Failed as e:
        print(f"failed: {_masked(str(e))}", file=sys.stderr)
        return EXIT_FAILED
    print("ok")
    print(side["output"]["sha256"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
