#!/usr/bin/env python
"""The barred branch's T3 file with the C2 factor of row R17 times the square of row K23's outside table: the
combination file CK of ledger row R18 (decision log 2026-11-21 (a)), fixed before any score of C2 or of the K23 file.

A thin builder, in the manner of build_two_knockout.py and build_outside_form.py, which it imports unchanged with the
frozen builder (build_submission.py) and outside_forms.py. It builds a barred-branch file exactly as
`build_submission.py build DECISIONS --branch barred --submission s1|s2` does, with one change, for this process only:
the factor the frozen loader returns for `transfer_mab21l2_all` is replaced by

  identity-B   build_two_knockout.combine('identity')              = phi_M ** x, the base file byte for byte
  gate-c2-B    build_two_knockout.combine('c2') * T ** 0.0         = C2 on B (build_two_knockout's c2-B) byte for byte
  gate-k23-B   build_two_knockout.combine('identity') * T ** 2.0   = the K23 file on B (build_outside_form.py
                                                                     xall-B --table-exponent 2.0) byte for byte
  ck-B         build_two_knockout.combine('c2') * T ** 2.0         = CK = phi_M ** (x c) * phi_G ** c * T ** 2

for B in s1, s1x1.5, s1x2.0, s2 (x = 1, 1.5, 2.0, 1). phi_M and phi_G are the two committed factor tables of the run
gata4_20261028T122425Z, loaded by the frozen loader and checked against build_two_knockout's pins; c is C2's own
c = SD(x ln phi_M) / SD(x ln phi_M + ln phi_G), computed by build_two_knockout.combine from those two tables alone and
never with T. T is the outside table `all`, loaded only through outside_forms.load_table (folder from VEC_T3_EXT_TABLES,
listing and table at the pins in outside_forms.py). The table exponent is TABLE_EXPONENT = 2.0, fixed here; T is not
multiplied by c. The product is computed in this order: part = combine(...)[0], then part * (T ** e). Nothing is
rescaled or clipped after it; t3_decompose.build receives it at w = 1 and applies it once.

Refused by name (exit 2): c1-*, c2-* (build_two_knockout.py's forms), xall-*, xunrel-* and K22's and K23's arm names
(build_outside_form.py's), and anything else. No number or path is taken from the command line: seed 1000001 and n
'max' are the only values accepted. Exit 2 also for: frozen files off their pins; a tree that is not clean; for a form
with T, pins not filled or the committed K23 reading not naming exactly x_e2 (build_outside_form.check_reading_k23);
for gate-c2-B and ck-B, experiments/harness/TRIES.md untracked or without row R17; for ck-B, without a row R18 that
names CK; the frozen decision-file check; a spec factor other than {mab21l2, 1.0}; a table off its pin; a folder,
listing or record that outside_forms refuses.

Exit 4 (ck-B only; the form's own refusal: CK is not built and nothing replaces it): |corr(ln phi_G, ln T)| over the
genes where T != 1 at or above outside_forms.MAX_ABS_CORR, or undefined; CK differing from the C2 factor at fewer than
outside_forms.MIN_CHANGED genes; SD(ln CK) > 1.05 * SD(x ln phi_M); CK not finite and positive.

The sidecar is the frozen builder's plus one block `combination_build`. It adds no cell number. `inputs.factor` stays
the frozen loader's record of phi_M; the factor applied is the one `combination_build` names (`inputs.factor_note`).

  VEC_T3_EXT_TABLES=~/v/vec-ext/t3/tables uv run --frozen --no-sync python -I -B \\
      experiments/t3_decomp/build_combination.py build DECISIONS.json --form ck-s1 --wt data/raw/T2/E8.75.h5ad \\
      --run-dir runs/t3_decomp/gata4_<STAMP> --board T3:bcatenin --seed 1000001 --out runs/<folder>/T3_ck_s1.h5ad
"""

from __future__ import annotations

import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import typer

import vec.provenance  # noqa: F401  (records this build's start, before the method code below is loaded)

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import build_outside_form as bo  # noqa: E402  (the K23 reading's check and the frozen check)
import build_submission as bs  # noqa: E402
import build_two_knockout as tk  # noqa: E402  (combine, the pins of phi_M and phi_G, the frozen and clean checks)
import outside_forms as ext  # noqa: E402  (the listing, the table, the figures)

SCRIPT = "experiments/t3_decomp/build_combination.py"
SCRIPT_PATH = Path(__file__).resolve()
BRANCH = "barred"
FROZEN_FACTOR = {"source": "mab21l2", "exponent": 1.0}
TABLE_M, TABLE_G, TABLE_T = tk.TABLE_M, tk.TABLE_G, "all"
TABLE_EXPONENT = 2.0  # K23's named exponent; fixed here, never from the command line
SD_BOUND = 1.05       # SD(ln CK) <= SD_BOUND * SD(x ln phi_M), K23's 5% bound
BASES = dict(tk.BASES)  # base -> (frozen slot, x, the reference sha256 of the base file built on the DGX)
# the parts' recorded hashes on each base (DGX builds): C2 from experiments/builds/results/t3two_fb01487.SHA256SUMS,
# the K23 file from experiments/builds/results/t3k23_49f3fd5.SHA256SUMS
C2_SHA256 = {
    "s1": "047506dbab19f007ce45820e3ab79c57b760fbbfc8da53058e4748b8c9930eac",
    "s1x1.5": "12a0dffbfc4fa0cdbedaef992edf17778e2610792008cc114e12b02569f8c1c8",
    "s1x2.0": "9559c9cd785b2eefc07847d4094d2d002daaac6a03bc2ba422eb00cd9e1d8d00",
    "s2": "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7",
}
K23_SHA256 = {
    "s1": "f3af008c12e6a3f29e20a3fc770a43671012923d8dc25fbfada29138f1d8854f",
    "s1x1.5": "0b31c03e301fca601957f06cbe1d31a204e425a287a70394e56b4a6ad9ee7934",
    "s1x2.0": "7e2815bd6f1d7c43b116e22e9976d0629823b15718dda5e1b3a9cde5b5ef82f5",
    "s2": "77ff0b947f7f409c9ce171bdc59a44de3d3c2ef2f3e5a7332011a459ce9ee9c3",
}
# form kind -> (build_two_knockout.combine's kind, the table exponent; None: no table)
KINDS = {"identity": ("identity", None), "gate-c2": ("c2", 0.0), "gate-k23": ("identity", TABLE_EXPONENT),
         "ck": ("c2", TABLE_EXPONENT)}
FORMS = tuple(f"{k}-{b}" for k in KINDS for b in BASES)
REFUSED_KINDS = {"c1": "a form of build_two_knockout.py (row R17), never built here",
                 "c2": "a form of build_two_knockout.py (row R17), never built here"}
REFUSED_KINDS |= {k: "a form or arm of build_outside_form.py (rows K22, K23), never built here"
                  for k in ("xall", "xunrel", *bo.REFUSED_KINDS)}
TRIES_PATH = HERE.parents[1] / "experiments" / "harness" / "TRIES.md"
ROW_R17, ROW_R18 = "| R17 (", "| R18 ("
WHERE = ("the frozen loader's phi_M and phi_G (each already rescaled and clipped to [0.25, 4] by gene_factor) give the "
         "C2 factor through build_two_knockout.combine; the outside table T (outside_forms.load_table: rescaled, "
         "clipped, exactly 1 at absent genes) enters as T ** e after it; t3_decompose.build receives the product at "
         "w = 1 and applies it through apply_gene_factor after the cells are drawn; nothing is rescaled or clipped "
         "after the product")
LIST_ENTRY = "docs/PLAN.md decision log 2026-11-21 (a); experiments/harness/TRIES.md row R18"
LIST_SEED, LIST_N_CELLS = 1_000_001, "max"
LIST_INPUTS_REFUSAL = "the list fixes seed 1000001 and max_cells"
FACTOR_NOTE = "the frozen loader's record of phi_M only; the factor applied is combination_build's"
EXIT_REFUSED, EXIT_FORM_REFUSED = 2, 4

Refused = bs.Refused


class FormRefused(Refused):
    """The form's own refusal (exit status 4): CK is not built and nothing replaces it."""


app = typer.Typer(no_args_is_help=True, add_completion=False)


def parse_form(form) -> tuple[str, str]:
    """--form: exactly one of FORMS; returns (kind, base), else Refused (the other builders' forms by name)."""
    f = str(form) if form is not None else ""
    for kind in REFUSED_KINDS:
        if f.startswith(f"{kind}-") or f == kind:
            raise Refused(f"--form {form!r} is refused: {kind} is {REFUSED_KINDS[kind]}")
    if f not in FORMS:
        raise Refused(f"--form {form!r} is not one of {list(FORMS)}")
    kind, base = f.rsplit("-", 1)
    return kind, base


def check_rows(kind: str, path: Path | None = None, check_tracked: bool = True) -> dict | None:
    """gate-c2 and ck: TRIES.md tracked and holding row R17; ck: also a row R18 that names CK. Returns a record."""
    if kind not in ("gate-c2", "ck"):
        return None
    path = TRIES_PATH if path is None else Path(path)
    if path.is_symlink() or not path.is_file():
        raise Refused(f"no ledger at {path.name}: not built")
    if check_tracked and not ext.tracked(path):
        raise Refused(f"{path.name} is not committed: not built")
    lines = path.read_text(encoding="utf-8").splitlines()
    if not any(ln.startswith(ROW_R17) for ln in lines):
        raise Refused("the ledger holds no row R17: not built")
    if kind == "ck" and not any(ln.startswith(ROW_R18) and "CK" in ln for ln in lines):
        raise Refused("the ledger holds no row R18 that names CK: CK is built only after its row")
    return {"file": "experiments/harness/TRIES.md", "rows": ["R17", "R18"] if kind == "ck" else ["R17"]}


def _sd(v: np.ndarray) -> float:
    return float(np.std(v, ddof=0))


def ck_factor(kind: str, x: float, phi_m, phi_g, T, e) -> tuple[np.ndarray, dict]:
    """tk.combine(kind, x, phi_m, phi_g)[0] * (T ** e), in that order (the part alone when T is None); returns the
    factor and a record. c is the part's own (combine's, from phi_M and phi_G alone)."""
    part, prec = tk.combine(kind, x, phi_m, phi_g)
    if T is None:
        if e is not None:
            raise Refused("a table exponent needs a table")
        return part, {"part": prec, "table_exponent": None}
    if e is None:
        raise Refused("a table needs its exponent")
    T = np.asarray(T, dtype=np.float64)
    if T.shape != part.shape:
        raise Refused("the outside table and the factor tables do not cover the same genes")
    phi = part * (T ** e)
    return phi, {"part": prec, "table_exponent": float(e)}


def own_refusals(x: float, phi_m, phi_g, T, part: np.ndarray, phi: np.ndarray) -> None:
    """ck's own refusals (FormRefused, exit 4); see the module."""
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise FormRefused("CK is not finite and positive")
    c = ext.table_figures(phi_g, T)["corr_ln_phi_m"]  # corr(ln phi_G, ln T) over the genes where T != 1
    if c is None or abs(c) >= ext.MAX_ABS_CORR:
        raise FormRefused(f"|corr(ln phi_G, ln T)| over T != 1 is {'undefined' if c is None else f'{abs(c):.3f}'} "
                          f"(>= {ext.MAX_ABS_CORR})")
    changed = int((phi != part).sum())
    if changed < ext.MIN_CHANGED:
        raise FormRefused(f"CK differs from the C2 factor at {changed} genes (< {ext.MIN_CHANGED})")
    sd_b, sd_ck = _sd(x * np.log(np.asarray(phi_m, dtype=np.float64))), _sd(np.log(phi))
    if sd_ck > SD_BOUND * sd_b:
        raise FormRefused(f"SD(ln CK) {sd_ck:.4f} > {SD_BOUND} * SD(x ln phi_M) {sd_b:.4f}")


def figures(x: float, phi_m, phi_g, T, part: np.ndarray, phi: np.ndarray) -> dict:
    """The sidecar's target-free figures of the product (none when there is no table)."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    sd_b, sd_c2, sd_ck = _sd(x * np.log(phi_m)), _sd(np.log(part)), _sd(np.log(phi))
    out = {"sd_ln_base": sd_b, "sd_ln_c2": sd_c2, "sd_ln_factor": sd_ck,
           "sd_ratio_factor_to_base": sd_ck / sd_b if sd_b > 0 else None, "sd_bound": SD_BOUND,
           "corr_ln_t_ln_phi_g": None, "corr_ln_t_ln_phi_m": None, "genes_changed_against_c2": None,
           "factor": None, "c2": None}
    if T is None or not (np.all(np.isfinite(phi) & (phi > 0))):
        return out
    out |= {"corr_ln_t_ln_phi_g": ext.table_figures(phi_g, T)["corr_ln_phi_m"],
            "corr_ln_t_ln_phi_m": ext.table_figures(phi_m, T)["corr_ln_phi_m"],
            "genes_changed_against_c2": int((phi != part).sum()),
            "factor": ext.figures(part, phi), "c2": ext.figures(part, part)}
    for k in ("factor", "c2"):
        out[k].pop("genes_changed")
    return out


@contextmanager
def combination_loader(kind: str, x: float, out: dict, env: dict | None = None):
    """bs.load_factor, for this process only: the frozen loader for phi_M and phi_G (each checked against
    build_two_knockout's pins too), the outside table through outside_forms for the forms that use it, then
    ck_factor. `out` receives the records."""
    original = bs.load_factor
    ckind, e = KINDS[kind]

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != tk.TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({tk.TABLE_SHA256[nm][:16]}): not built")
        T, trec = None, None
        if e is not None:
            try:
                listing, lsha = ext.read_listing()
                T, trec = ext.load_table(TABLE_T, list(genes), listing, lsha, env=env)
            except ext.TableRefused as err:
                raise Refused(str(err)) from None
        phi, frec = ck_factor(ckind, x, phi_m, phi_g, T, e)
        part = tk.combine(ckind, x, phi_m, phi_g)[0]
        if kind == "ck":
            own_refusals(x, phi_m, phi_g, T, part, phi)
        out.update(table_m=rec_m, table_g=rec_g, table=trec, factor=frec,
                   figures=figures(x, phi_m, phi_g, T, part, phi), n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


def build_form(decision: Path, form, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None = None,
               seed: int = 1_000_001, overwrite: bool = False, threads: int = 16, env: dict | None = None,
               tries: Path | None = None) -> dict:
    """The named form through bs.build_submission; returns the sidecar (the frozen builder's plus
    `combination_build`)."""
    kind, base = parse_form(form)
    slot, x, base_sha = BASES[base]
    ckind, e = KINDS[kind]
    frozen = tk.check_frozen()
    bo.check_frozen()
    tree = tk.check_clean()
    k23rec = None
    if e is not None:
        try:
            ext.check_pins(ext.K22_LISTING_SHA256, ext.K22_TABLE_SHA256)
        except ext.TableRefused as err:
            raise Refused(str(err)) from None
        k23rec = bo.check_reading_k23(TABLE_EXPONENT)
    rows = check_rows(kind, tries)
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    got: dict = {}
    with combination_loader(kind, x, got, env):
        side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                   run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the form's factor once at w = 1")
    side["inputs"]["factor_note"] = FACTOR_NOTE
    prec, trec = got["factor"]["part"], got["table"]
    tables = {k: {"file": v["file"], "sha256": v["sha256"], "recorded_in": v["recorded_in"],
                  "pinned_sha256": tk.TABLE_SHA256[k]} for k, v in ((TABLE_M, got["table_m"]),
                                                                    (TABLE_G, got["table_g"]))}
    tables[TABLE_T] = None if trec is None else {
        k: trec[k] for k in ("name", "file", "dir_realpath", "sha256", "listing_sha256", "pinned_sha256",
                             "combine_sha256", "manifest_sha256", "studies", "tiers", "sign_rule", "refstats_sha256")}
    side["combination_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH), "form": f"{kind}-{base}",
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha,
                 "rng_streams": side["inputs"]["rng_streams"]},
        "x": x, "c": prec["c"], "gata4_weight": prec["gata4_weight"], "mab21l2_exponent": prec["mab21l2_exponent"],
        "table_exponent": e, "part": prec, "tables": tables,
        "listing_sha256": None if trec is None else trec["listing_sha256"],
        "listing_pinned_sha256": None if trec is None else ext.K22_LISTING_SHA256,
        "k23_reading": k23rec, "ledger_rows": rows,
        "parts_sha256": {"c2": C2_SHA256[base], "k23": K23_SHA256[base]},
        "figures": got["figures"], "where": WHERE,
        "frozen_files_sha256": frozen, "tree_at_start": tree, "list_entry": LIST_ENTRY,
        "identity": {"identity": f"the {base} base file byte for byte ({base_sha[:12]})",
                     "gate-c2": f"C2 on {base} byte for byte ({C2_SHA256[base][:12]})",
                     "gate-k23": f"the K23 file on {base} byte for byte ({K23_SHA256[base][:12]})"}.get(kind),
        "depth_step": None,
    }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


@app.callback()
def main() -> None:
    """The barred file with C2's factor times the square of the outside table, in one of the fixed forms."""


@app.command()
def build(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    form: str = typer.Option(..., help=f"one of {', '.join(FORMS)}; nothing else"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the frozen builder"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    out: Path = typer.Option(..., help="output .h5ad (sidecar: same name with .json)"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(1_000_001, help="1000001, the list's seed; any other is refused"),
    n_cells: str = typer.Option("max", help="'max', the list's count: the board's max_cells, read in code, never "
                                "printed; any other is refused"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    overwrite: bool = typer.Option(False, help="replace an existing output"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """Build one fixed form."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        parse_form(form)
        b = bs.board_spec(board, boards_dir)
        s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                       overwrite=overwrite, threads=threads)
    except FormRefused as e:
        bs.say(f"refused (the form's own refusal): {e}", "red")
        raise typer.Exit(EXIT_FORM_REFUSED) from None
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    chk = s["output"]["check"]
    cb = s["combination_build"]
    bs.say(f"form {cb['form']}: base {cb['base']['name']} (x = {cb['x']}), c = {cb['c']:.4f}, table exponent "
           f"{cb['table_exponent']}")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


if __name__ == "__main__":
    app()
