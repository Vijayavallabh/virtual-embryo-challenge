#!/usr/bin/env python
"""Row R72's sign cube files: C2 (`build_two_knockout.py --form c2-s2`, recorded sha256 47111ac1d07e...) with one
change, its factor replaced by a vertex of the cube (sign_cube.py; T3_MOON5_design.md 2.2 with fixes F1 to F10).

A thin builder, in build_two_knockout.py's manner. It imports the frozen builder (build_submission.py) and script
(t3_decompose.py) unchanged and builds the barred s2 exactly as C2's command does (its spec, composition, seed
1000001, max_cells read by code, streams "submission:barred:s2" and "...:assemble", coordinates rule), with one
change: the factor the frozen loader returns is replaced, for this process only, by the vertex's factor

  ln phi_v = ln phi_C2 + sum_j (k_v s_j - 1) Q(X_j),   j = S (stage), M (Mab21l2), G (Gata4)

on C2's factor phi_C2 = combine('c2', 1.0) of the two pinned tables, with the axes of the four committed WT-only tables
(`cube_tables.py tables`; MANIFEST.json) and Q projecting off B's columns m_L * phi_C2 for every label with a cell in
C2's draw (the rows C2's own draw takes, counts as drawn, computed here on a fresh generator before the build; row
R72's every label); then sign_cube.flip_factor's SD hold, clip, den constant (WT E8.75's mean within-cell proportions,
build_outside_form.wt_den) and Newton hold. The board's own writer writes the file (bs.build_submission:
vec.io.write_submission), and its stamp follows.

  uv run --frozen --no-sync python -I -B experiments/t3_decomp/build_sign_cube.py build DECISIONS.json \\
      --form c2-s2 --vertex <v> --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> \\
      --board T3:bcatenin --seed 1000001 --cube-tables experiments/t3_decomp/results/t3_cube/MANIFEST.json \\
      --captain-ruling R72 --out runs/<folder>/T3_cube_<v>.h5ad

The seven files of row R72 are the vertices ppn, pnp, pnn, npp, npn, nnp, nnn. `ppp` is the identity gate: its factor
is C2's own array (every coefficient 0), so it must give C2's recorded bytes; it is never a file to upload.
Refusals (exit 2, before anything is built): a form other than c2-s2 (taken as typed); a vertex other than the eight;
a seed other than 1000001 or an --n-cells other than 'max'; no --captain-ruling, or any ruling but R72; a ledger
(experiments/harness/TRIES.md) without exactly one row R72 whose first cell is the captain's order ("ruling, no run",
build_outside_form.captain_order_names) naming `--captain-ruling R72` and, for each of the seven, the vertex's option
string `--form c2-s2 --vertex <v>` as one whole backtick-delimited string; a manifest whose sha256 is not the one row
K31 records (k31_cube_read.row_k31), or tables whose bytes are not the manifest's; sign_cube.py or cube_tables.py
not the bytes the manifest records (the code that made the tables, ran G0 and row K31's read); frozen files, pinned
tables or the decision file off their pins; a tree that is not clean. Exit 6 (nothing written): the vertex's factor
cannot be made (k_v outside [0.5, 2] or no root, the hold not converged) or a held label's expected library is off
C2's by more than 1% (row R72's words). Exit 5 (nothing written): the factor the build would apply is not the one
checked.
  vertices --reading R --g0 G   the vertices row K31's reading says to build, after checking it against G0's
                                reading (its sha256, CONTINUE, every vertex one G0 leaves): prints them, exit 0;
                                exit 3 when the reading builds nothing; exit 2 refused
"""

from __future__ import annotations

import hashlib
import json
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import typer

import vec.provenance  # noqa: F401  (records this build's start, before the method code below is loaded)

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import build_submission as bs  # noqa: E402
import build_two_knockout as btk  # noqa: E402
import cube_tables as ct  # noqa: E402
import sign_cube as sc  # noqa: E402

SCRIPT = "experiments/t3_decomp/build_sign_cube.py"
SCRIPT_PATH = Path(__file__).resolve()
FORM = "c2-s2"
ROW = "R72"
K31_ROW = "K31"
BRANCH, SLOT, STREAM = "barred", "s2", "submission:barred:s2"
SEED, N_CELLS = 1_000_001, "max"
TRIES_FILE = btk.TRIES_FILE
C2_S2_SHA256 = btk.C2_S2_SHA256  # C2's file, t3two_fb01487.SHA256SUMS
FILES = {v: f"T3_cube_{v}.h5ad" for v in sc.NAMED}
EXIT_REFUSED, EXIT_FACTOR_MISMATCH, EXIT_GUARD, EXIT_NOTHING = 2, 5, 6, 3
LINE = ("row R72's sign cube on C2: ln phi_v = ln phi_C2 + sum_j (k_v s_j - 1) Q(X_j), the axes "
        "X_S = c (S_M - Y_M) / 2, X_M = c (l_M - S_M), X_G = c (l_G - (S_G + Y_G) / 2) from the committed WT-only "
        "tables, Q off B's columns m_L * phi_C2 (every label with a cell in C2's draw), k_v for SD(ln phi_C2), the "
        "clip to [min(0.25, phi_C2), max(4, phi_C2)], the den constant and the Newton hold of every label's expected "
        "library, each within 1% of C2's")

Refused = bs.Refused


class FactorMismatch(Refused):
    """The factor the build applies is not the one the checks describe (exit 5)."""


class GuardFailed(Refused):
    """The vertex's factor cannot be made or breaks the library band (exit 6, nothing written)."""


app = typer.Typer(no_args_is_help=True, add_completion=False)


def option_string(vertex: str) -> str:
    return f"--form {FORM} --vertex {vertex}"


def check_command(form, vertex, captain_ruling, seed=SEED, n_cells=N_CELLS, tries: Path | None = None) -> dict:
    """The command's refusals that need no data (Refused, exit 2): the form, the vertex, the seed and count, the
    ruling, and row R72's first cell naming the vertex's whole option string (ppp, the identity gate, needs the order
    only)."""
    if form != FORM:
        raise Refused(f"--form {form!r} is refused: only {FORM}")
    if not isinstance(vertex, str) or vertex not in sc.VERTICES:
        raise Refused(f"--vertex {vertex!r} is refused: one of {list(sc.VERTICES)}")
    if seed != SEED or n_cells != N_CELLS:
        raise Refused("C2's command fixes seed 1000001 and max_cells")
    if captain_ruling != ROW:
        raise Refused(f"--captain-ruling {captain_ruling!r} is refused: only row {ROW} orders the cube's files")
    p = TRIES_FILE if tries is None else Path(tries)
    first = btk.ledger_first_cell(ROW, p)
    if not btk._outside_form().captain_order_names(ROW, f"`--captain-ruling {ROW}`", p) or first is None:
        raise Refused(f"experiments/harness/TRIES.md has no single row {ROW} whose first cell is the captain's order "
                      f"(\"ruling, no run\") naming `--captain-ruling {ROW}`")
    if vertex != sc.IDENTITY and not btk.names_whole(first, option_string(vertex)):
        raise Refused(f"row {ROW}'s first cell does not name `{option_string(vertex)}` as one whole string")
    return {"row": ROW, "vertex": vertex, "option": f"`{option_string(vertex)}`",
            "file": FILES.get(vertex, "identity gate (never uploaded)"), "ledger": "experiments/harness/TRIES.md",
            "ledger_sha256": bs.sha256(p)}


def check_manifest(manifest: Path, tries: Path | None = None) -> str:
    """The manifest's sha256 must be the one row K31 records; returns it (Refused otherwise)."""
    import k31_cube_read as kr

    want = kr.row_k31(TRIES_FILE if tries is None else Path(tries))
    got = bs.sha256(Path(manifest)) if Path(manifest).is_file() else None
    if got != want:
        raise Refused(f"the cube tables manifest's sha256 is not the one row {K31_ROW} records")
    return got


def check_code(man: dict) -> dict:
    """The cube's algebra and tables script must be the bytes that made the tables (the manifest's record), so that
    the files are made by the code G0 and row K31's read ran with (Refused otherwise)."""
    got = {"sign_cube.py": bs.sha256(HERE / "sign_cube.py"), "cube_tables.py": bs.sha256(HERE / "cube_tables.py")}
    want = {"sign_cube.py": (man.get("code_sha256") or {}).get("sign_cube.py"),
            "cube_tables.py": man.get("script_sha256")}
    if got != want:
        raise Refused("sign_cube.py or cube_tables.py is not the code that made the cube tables (the manifest's "
                      "record): not built")
    return got


def drawn_rows(ref, kw: dict, n_cells: int, seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """The rows C2's own draw takes (construct with C2's keywords and no factor, on a fresh stable_rng(seed, STREAM):
    the rows are drawn before the factor step) and their labels."""
    built = bs.construct(ref, n_cells, bs.t3d.stable_rng(seed, STREAM), **{**kw, "phi": None})
    if not np.all(built.source == "wt"):
        raise Refused("C2's cells come from the reference pool only")
    return built.row, built.lab


def vertex_of(phi_m, phi_g, tabs: dict, counts, lab, den, vertex: str) -> tuple[np.ndarray, dict, np.ndarray]:
    """The vertex's factor from the two pinned tables, the four cube tables (phi, gene order of WT E8.75), C2's drawn
    counts and labels, and den: (phi_v, record, phi_C2). Raises sign_cube.CubeRefused when it cannot be made."""
    phi_c2, crec = btk.combine("c2", 1.0, phi_m, phi_g)
    c = float(crec["c"])
    lt = {t: np.log(np.asarray(tabs[t], dtype=np.float64)) for t in ct.TABLES}
    ax = sc.cube_axes(np.log(phi_m), np.log(phi_g), lt["S_M"], lt["Y_M"], lt["S_G"], lt["Y_G"], c)
    held, M = sc.label_rows(counts, lab)
    B = sc.basis(M, phi_c2)
    Q = sc.projector(B)
    q = {j: Q(ax[j]) for j in sc.AXES}
    phi, rec = sc.vertex_factor(phi_c2, q, vertex, M, den, B)
    rec |= {"c": c, "held_labels": int(len(held)), "sd_q": {j: sc.sd(q[j]) for j in sc.AXES},
            "c2_factor_sha256": sc.factor_sha256(phi_c2), "identity": bool(np.array_equal(phi, phi_c2))}
    return phi, rec, phi_c2


@contextmanager
def vertex_loader(phi: np.ndarray, genes: list[str], out: dict):
    """bs.load_factor, for this process only: the frozen loader for phi_M (checked against the pinned table), then the
    vertex's factor in its place. FactorMismatch unless the build's genes are those the factor was made over."""
    original = bs.load_factor
    sha = sc.factor_sha256(phi)

    def load(run_dir, name, g, dec):
        if name != btk.TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {btk.TABLE_M!r}: not built")
        if list(g) != list(genes):
            raise FactorMismatch("the build's genes are not those the vertex's factor was made over: not built")
        _, rec_m = original(run_dir, btk.TABLE_M, g, dec)
        if rec_m["sha256"] != btk.TABLE_SHA256[btk.TABLE_M]:
            raise Refused(f"{rec_m['file']}: not the pinned table: not built")
        if sc.factor_sha256(phi) != sha:
            raise FactorMismatch("the factor changed after its checks: not built")
        out.update(n_loads=out.get("n_loads", 0) + 1, factor_sha256=sha, table_m=rec_m)
        return np.array(phi, dtype=np.float64, copy=True), rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


def build_vertex(decision: Path, form, vertex, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None,
                 manifest: Path, captain_ruling=None, seed: int = SEED, overwrite: bool = False, threads: int = 16,
                 tries: Path | None = None) -> dict:
    """One vertex through bs.build_submission; returns the sidecar (written next to the file)."""
    order = check_command(form, vertex, captain_ruling, seed, N_CELLS, tries)
    man_sha = check_manifest(Path(manifest), tries)
    frozen = btk.check_frozen()
    tree = btk.check_clean()
    out_p = Path(out)
    if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
        raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    spec = bs.barred_table(dec, rules)[SLOT]
    if spec["factor"] != btk.FROZEN_FACTOR:
        raise Refused(f"the frozen barred s2 factor is {spec['factor']}, not {btk.FROZEN_FACTOR}: not built")
    genes, den = btk._outside_form().wt_den(Path(wt))
    genes, den = list(genes), np.asarray(den, dtype=np.float64)
    tabs, man, _ = ct.load_tables(Path(manifest), genes, man_sha)
    code = check_code(man)
    phi_m, rec_m = bs.load_factor(run_dir, btk.TABLE_M, genes, dec)
    phi_g, rec_g = bs.load_factor(run_dir, btk.TABLE_G, genes, dec)
    for nm, rec in ((btk.TABLE_M, rec_m), (btk.TABLE_G, rec_g)):
        if rec["sha256"] != btk.TABLE_SHA256[nm]:
            raise Refused(f"{rec['file']}: not the pinned table: not built")
    got: dict = {}
    with btk.combined_loader("c2", 1.0, got):
        pf = bs.preflight(Path(decision), BRANCH, SLOT, board, n_cells)
        ref = bs.load_reference(Path(wt), pf["preset"], rules, None, threads)
        kw, _ = bs.inputs_for(pf["spec"], pf["dec"], ref, run_dir)
    if list(ref.genes) != list(genes):
        raise Refused("the reference's genes are not WT E8.75's")
    rows, lab = drawn_rows(ref, kw, n_cells, seed)
    try:
        phi, vrec, phi_c2 = vertex_of(phi_m, phi_g, tabs, ref.pool.C[rows], lab, den, vertex)
    except sc.CubeRefused as e:
        raise GuardFailed(f"vertex {vertex}: {e}: not built") from None
    if not np.array_equal(kw["phi"], phi_c2):
        raise FactorMismatch("C2's command does not apply C2's factor: not built")
    if not sc.within_band(vrec):
        raise GuardFailed(f"vertex {vertex}: a held label's expected library is off C2's by "
                          f"{vrec['max_abs_library_dev']:.3g}, beyond {sc.LIB_BAND:g}: not built")
    hook: dict = {}
    draws: list = []
    with vertex_loader(phi, genes, hook), bs.observe_draws(draws):
        side = bs.build_submission(Path(decision), BRANCH, SLOT, Path(wt), board, n_cells, out_p, run_dir=run_dir,
                                   seed=seed, overwrite=overwrite, threads=threads)
    if hook.get("n_loads") != 1 or hook.get("factor_sha256") != vrec["factor_sha256"]:
        raise RuntimeError("the build did not apply the vertex's checked factor exactly once")
    built_rows = np.concatenate([np.asarray(r, dtype=np.int64) for _, r in draws if len(r)])
    if not np.array_equal(built_rows, np.asarray(rows, dtype=np.int64)):
        raise RuntimeError("the build did not draw C2's rows (the rows B was made from)")
    if side["construction"]["factor_exponent"] != 1.0 or side["spec"]["factor"] != btk.FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the vertex's factor at w = 1 on C2's spec")
    side["inputs"]["factor_note"] = "the frozen loader's record of phi_M only; the factor applied is cube_build's"
    side["cube_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH), "form": FORM, "vertex": vertex,
        "captain_ruling": order, "line": LINE, "manifest": str(manifest), "manifest_sha256": man_sha,
        "tables_sha256": {t: r["sha256"] for t, r in man["tables"].items()},
        "pinned_tables": {btk.TABLE_M: rec_m["sha256"], btk.TABLE_G: rec_g["sha256"]},
        "factor": {k: v for k, v in vrec.items() if k not in ("k_solve",)} | {"k_solve": vrec["k_solve"]},
        "factor_sha256_applied": hook["factor_sha256"], "c2_reference_sha256": C2_S2_SHA256,
        "identity_gate": vertex == sc.IDENTITY, "upload": vertex != sc.IDENTITY,
        "unchanged": "C2's spec, composition, rows (drawn before the factor step; checked equal to the rows B was "
                     "made from), count (max_cells, read by code), seed, streams and coordinates rule",
        "frozen_files_sha256": frozen, "tree_at_start": tree, "cube_code_sha256": code,
    }
    bs.dump_json(side, out_p.with_suffix(".json"))
    return side


def vertices_to_build(reading: Path, g0: Path) -> list[str]:
    """Row K31's reading's vertex list, checked against G0's reading: the reading names g0.json's sha256, both say
    CONTINUE, and every vertex is one of the seven that G0 leaves (Refused otherwise)."""
    r = json.loads(Path(reading).read_text())
    raw = Path(g0).read_bytes()
    g = json.loads(raw)
    if (r.get("g0") or {}).get("sha256") != hashlib.sha256(raw).hexdigest():
        raise Refused("the reading was not made from this G0 reading")
    if r.get("outcome") != "CONTINUE":
        return []
    if g.get("outcome") != "CONTINUE":
        raise Refused("the reading says CONTINUE but G0 did not")
    out = list(r.get("build") or [])
    if any(v not in sc.NAMED or v not in (g.get("vertices_left") or []) for v in out) or len(set(out)) != len(out):
        raise Refused("the reading names a vertex G0 does not leave")
    return out


@app.callback()
def main() -> None:
    """Row R72's sign cube on C2: one vertex per call."""


@app.command()
def build(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    form: str = typer.Option(..., help="c2-s2 only"),
    vertex: str = typer.Option(..., help="ppp (the identity gate) or one of the seven: ppn pnp pnn npp npn nnp nnn"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the frozen builder"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    out: Path = typer.Option(..., help="output .h5ad (sidecar: same name with .json)"),
    cube_tables: Path = typer.Option(..., help="the cube tables' MANIFEST.json (its sha256: row K31's)"),
    captain_ruling: str = typer.Option(None, help="R72 only"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(SEED, help="1000001 only"),
    n_cells: str = typer.Option(N_CELLS, help="'max' only: the board's max_cells, read in code, never printed"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    overwrite: bool = typer.Option(False, help="replace an existing output"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """Build one vertex of the cube."""
    try:
        check_command(form, vertex, captain_ruling, seed, n_cells)
        check_manifest(cube_tables)
        b = bs.board_spec(board, boards_dir)
        s = build_vertex(decision, form, vertex, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir, cube_tables,
                         captain_ruling, seed, overwrite, threads)
    except GuardFailed as e:
        bs.say(f"refused (the vertex's checks, before anything is written): {e}", "red")
        raise typer.Exit(EXIT_GUARD) from None
    except FactorMismatch as e:
        bs.say(f"refused (the factor is not the one checked): {e}", "red")
        raise typer.Exit(EXIT_FACTOR_MISMATCH) from None
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    chk = s["output"]["check"]
    f = s["cube_build"]["factor"]
    bs.say(f"vertex {vertex}: k {f['k']:.4f}, hold rounds {f['hold']['rounds']}, largest label library "
           f"{f['max_abs_library_dev']:.2e}" + (" (the identity gate)" if vertex == sc.IDENTITY else ""))
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


@app.command()
def vertices(
    reading: Path = typer.Option(..., help="row K31's reading.json"),
    g0: Path = typer.Option(..., help="G0's g0.json"),
):
    """Prints the vertices row K31's reading builds (space-separated); exit 3 when none, 2 refused."""
    try:
        vs = vertices_to_build(reading, g0)
    except (Refused, OSError, ValueError) as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    if not vs:
        raise typer.Exit(EXIT_NOTHING)
    print(" ".join(vs))
    raise typer.Exit(0)


if __name__ == "__main__":
    app()
