#!/usr/bin/env python
"""The barred branch's s1 with its Mab21l2 factor raised to a named exponent (docs/PLAN.md decision log 2026-11-20 (h)).

A thin builder. It imports the frozen builder (build_submission.py) and script (t3_decompose.py) unchanged and builds
the barred s1 exactly as `build_submission.py build DECISIONS --branch barred --submission s1` does, with one change:
the s1 factor's exponent in the spec that rules.yaml note 3 gives (`build_submission.barred_table`, exponent 1) is
replaced by --exponent, for this process only. Everything else is the frozen builder's own code: the decision file and
its hash check, the central composition at m = D7m, the drawn cells and their stream ("submission:barred:s1"), the
row order and coordinates (stream "...:assemble"), the factor file and its sha256 check, no depth, stage or coverage
step, the seed rule, the board check, the file, its build stamp and its sidecar.

Where the exponent enters. `t3_decompose.build` applies `phi ** w` (w = the spec's exponent) through
`apply_gene_factor`, after the cells are drawn and before lognorm. phi is the factor file
`<run>/factors/transfer_mab21l2_all.csv` as loaded, whose values `t3_decompose.gene_factor` already clipped to
[0.25, 4] when `transfer_factors` wrote them. So the exponent acts on the clipped factor and nothing clips again: at
exponent 2 the values lie in [0.0625, 16]. This is how the practice script of try K20 sets its arms (branch
campaign/t3-k20-exponent, `k20_barred_exponent.build_arm`: the barred keywords with `w` replaced) and how the frozen
builder itself applies any exponent (the allowed branch's 0.5 of rules.yaml notes 6 and 8). At w = 1 build passes phi
unchanged (no power is taken), so exponent 1.0 gives the frozen barred s1 byte for byte.

Only the exponents 1.0, 1.5 and 2.0 are accepted (1.0 is the identity check; 1.5 and 2.0 are the candidates fixed in
the decision log before any score of theirs existed). No value here is computed from a test score.

Refusals (exit status 2), before anything loads: another exponent; t3_decompose.py, rules.yaml or build_submission.py
whose sha256 differs from the bytes pinned below (main at ca2ebc5, the same bytes as at the third rehearsal,
dc7b6ae); a tree that is not clean by vec.provenance's rule (the rule the build stamp records: no tracked change and
no untracked file outside results/, runs/ and experiments/**/results/; build_all.sh refuses tracked changes); then
every refusal of the frozen builder, with --allow-hash-mismatch and --allow-unregistered-protocol never passed (so its
own check of the decision file's script and rules hashes against the files on disk must pass).

The sidecar is the frozen builder's, plus one block `exponent_build` (this script, its sha256, the exponent, the
frozen spec's factor, where the exponent enters, the pinned hashes). It adds no cell number.

  uv run --frozen --no-sync python -I -B experiments/t3_decomp/build_barred_exponent.py build DECISIONS.json \\
      --exponent 1.5 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> --board T3:bcatenin \\
      --seed 1000001 --out runs/<folder>/barred_s1_x1.5.h5ad
"""

from __future__ import annotations

import copy
import sys
from contextlib import contextmanager
from pathlib import Path

import typer

import vec.provenance  # noqa: F401  (records this build's start, before the method code below is loaded)

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import build_submission as bs  # noqa: E402

from vec import provenance  # noqa: E402

SCRIPT = "experiments/t3_decomp/build_barred_exponent.py"
SCRIPT_PATH = Path(__file__).resolve()
BRANCH, SUBMISSION = "barred", "s1"
EXPONENTS = (1.0, 1.5, 2.0)
FROZEN_FACTOR = {"source": "mab21l2", "exponent": 1.0}
# the frozen files' bytes on main at ca2ebc5, unchanged since the third rehearsal (dc7b6ae) that built the reference
# barred s1 (sha256 3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2)
FROZEN_SHA256 = {
    "t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
    "rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13",
    "build_submission.py": "1194a63993a27dde4dedfdb9180de035c5f8f2a57aaa7d97b38e7a15ea9e0c74",
}
WHERE = ("t3_decompose.build applies phi ** w (w = this exponent) through apply_gene_factor after the cells are drawn; "
         "phi is the factor file as loaded, already clipped to [0.25, 4] by gene_factor when transfer_factors wrote "
         "it; no clip after the power (as try K20's practice arms and the frozen builder's own exponents)")
CANDIDATES = "docs/PLAN.md decision log 2026-11-20 (h)"

Refused = bs.Refused
app = typer.Typer(no_args_is_help=True, add_completion=False)


def frozen_paths() -> dict[str, Path]:
    """The three frozen files as the frozen builder itself names them."""
    return {"t3_decompose.py": bs.T3D_PATH, "rules.yaml": bs.RULES_PATH, "build_submission.py": bs.BUILDER_PATH}


def parse_exponent(x) -> float:
    """--exponent: exactly one of EXPONENTS, else Refused."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        raise Refused(f"--exponent {x!r} is not a number; accepted: {list(EXPONENTS)}") from None
    if v not in EXPONENTS:
        raise Refused(f"--exponent {x!r} is not one of {list(EXPONENTS)} (decision log 2026-11-20 (h): 1.5 and 2.0; "
                      "1.0 is the identity check)")
    return v


def check_frozen(paths: dict[str, Path] | None = None) -> dict[str, str]:
    """The frozen files' sha256 on disk against the pinned bytes; returns them."""
    paths = frozen_paths() if paths is None else paths
    got = {k: bs.sha256(p) for k, p in paths.items()}
    bad = [k for k in FROZEN_SHA256 if got.get(k) != FROZEN_SHA256[k]]
    if bad:
        raise Refused(f"frozen file(s) {bad} differ from the pinned bytes (main at ca2ebc5): not built")
    return got


def check_clean() -> dict:
    """The tree by vec.provenance's rule (the one the build stamp records); refuses unless clean."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    if commit == "unknown" or unclean is None:
        raise Refused("git cannot tell the tree's state: not built (the stamp would not say clean)")
    if unclean:
        raise Refused(f"the tree is not clean ({len(unclean)} path(s), e.g. {unclean[:3]}): commit first (the build "
                      "stamp must say clean)")
    return {"commit": commit, "clean": True}


@contextmanager
def barred_s1_at(x: float):
    """bs.barred_table, for this process only, with s1's factor exponent set to x. Refuses unless the frozen table's
    s1 factor is the Mab21l2 factor at exponent 1; s2 and every other field are left as the frozen code gives them."""
    original = bs.barred_table

    def at_x(dec, rules):
        tab = original(dec, rules)
        if tab["s1"]["factor"] != FROZEN_FACTOR:
            raise Refused(f"the frozen barred s1 factor is {tab['s1']['factor']}, not {FROZEN_FACTOR}: not built")
        tab = copy.deepcopy(tab)
        tab["s1"]["factor"]["exponent"] = float(x)
        return tab

    bs.barred_table = at_x
    try:
        yield
    finally:
        bs.barred_table = original


def build_exponent(decision: Path, exponent, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None = None,
                   seed: int = 1_000_001, overwrite: bool = False, threads: int = 16) -> dict:
    """The barred s1 at the named exponent through bs.build_submission; returns the sidecar (written next to the
    file, the frozen builder's plus `exponent_build`)."""
    x = parse_exponent(exponent)
    frozen = check_frozen()
    tree = check_clean()
    rules = bs.load_rules()
    bs.check_hashes(bs.load_decision(Path(decision), rules), allow_mismatch=False)  # the frozen builder's own check
    with barred_s1_at(x):
        side = bs.build_submission(Path(decision), BRANCH, SUBMISSION, Path(wt), board, n_cells, Path(out),
                                   run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    if side["spec"]["factor"] != {"source": "mab21l2", "exponent": x} or side["construction"]["factor_exponent"] != x:
        raise RuntimeError("the built spec does not carry the named exponent")
    side["exponent_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH), "exponent": x,
        "frozen_spec_factor": FROZEN_FACTOR, "where": WHERE, "candidates": CANDIDATES,
        "frozen_files_sha256": frozen, "tree_at_start": tree,
        "identity": "exponent 1.0 is the frozen barred s1 byte for byte" if x == 1.0 else None,
    }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


@app.callback()
def main() -> None:
    """The barred s1 with its Mab21l2 factor at exponent 1.0, 1.5 or 2.0."""


@app.command()
def build(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    exponent: str = typer.Option(..., help="1.0 (identity check), 1.5 or 2.0; nothing else"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the frozen builder"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    out: Path = typer.Option(..., help="output .h5ad (sidecar: same name with .json)"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(1_000_001, help="as the frozen builder (never below vec.harness.BUILD_SEED_OFFSET)"),
    n_cells: str = typer.Option("max", help="'max' (default): the board's max_cells, read in code, never printed"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    overwrite: bool = typer.Option(False, help="replace an existing output"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """Build the barred s1 at the named exponent."""
    try:
        parse_exponent(exponent)
        b = bs.board_spec(board, boards_dir)
        s = build_exponent(decision, exponent, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir,
                           seed=seed, overwrite=overwrite, threads=threads)
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    chk = s["output"]["check"]
    comp = s["spec"]["composition"]
    bs.say(f"barred s1, factor mab21l2 at exponent {s['exponent_build']['exponent']}; composition {comp['kind']}")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


if __name__ == "__main__":
    app()

