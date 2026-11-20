#!/usr/bin/env python
"""The barred branch's T3 file with its Mab21l2 factor multiplied by one of row K24's scaled outside-knockout tables:
the one form K24's committed reading may name, on the base B that decision log 2026-11-20 (q) fixes. Nothing here
chooses a form, a table, an exponent or a base.

A thin builder in the manner of build_outside_form.py, which it imports and does not change (tests pin its bytes for
the CK build). It imports the frozen builder (build_submission.py) unchanged and builds a barred-branch file exactly as
`build_submission.py build DECISIONS --branch barred --submission s1|s2` does, with one change, for this process only:
the factor the frozen loader returns for `transfer_mab21l2_all` (phi_M, sha256-checked by the frozen loader against
the decision file and here against the pin) is replaced by

  k24sc3-B     phi_M ** x * sc3 ** e_sc3
  k24bulk-B    phi_M ** x * bulk ** e_bulk
  k24joint-B   phi_M ** x * sc3 ** (c e_sc3) * bulk ** (c e_bulk)

for B in s1, s1x1.5, s1x2.0, s2 (x = 1, 1.5, 2.0, 1), through k24_wider_tables.k24_table, the function K24's arms
apply. The tables are loaded by k24_wider_tables.load_k24_table: by key, from the folder VEC_T3_EXT_TABLES names, at the
sha256 of the committed experiments/t3_ext/TABLES_K24.sha256 and of the pins in k24_wider_tables.py. e_T and c are read
from the committed K24 reading's pre-gate (full precision) and recomputed from the loaded tables to a relative 1e-12.
Nothing is rescaled or clipped after the product.

Refused (exit 2, nothing built): any other form (the controls, the twins, x_e2, barred_form and identity by name);
another seed or n; frozen files off their pins; a tree that is not clean; pins not filled or the listing untracked or
off its pin; the folder variable unset or refused; bytes or records that disagree; a committed K24 reading
(experiments/t3_decomp/results/k24_wider_tables/reading.json, tracked by git) that is missing, VOID or incomplete, does
not name exactly the form's arm, records other tables than the listing and pins, holds the arm dropped, or whose e_T
or c differ from the ones recomputed. Exit 4: the form's own refusal (the product cannot be made, or it changes fewer
than 10 genes against phi_M ** x); nothing replaces it.

The sidecar is the frozen builder's plus one block `outside_build`: this script and its sha256, the form, the base, x,
each table's record with its register accessions (EXTERNAL_DATA.md rows found by accession), e_T, c, the reading's
sha256, priced_at 'x = 1 only' and the list entry.

  VEC_T3_EXT_TABLES=~/v/vec-ext/t3/tables uv run --frozen --no-sync python -I -B \\
      experiments/t3_decomp/build_k24_form.py build DECISIONS.json --form k24sc3-s1 --wt data/raw/T2/E8.75.h5ad \\
      --run-dir runs/t3_decomp/gata4_<STAMP> --board T3:bcatenin --seed 1000001 --out runs/<folder>/T3_k24sc3_s1.h5ad
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import typer

import vec.provenance  # noqa: F401  (records this build's start, before the method code below is loaded)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import build_outside_form as bo  # noqa: E402  (unchanged: its checks, bases and the frozen builder)
import k24_wider_tables as k24  # noqa: E402  (the tables, the size rule and the scaled table)
import outside_forms as ext  # noqa: E402

bs = bo.bs
SCRIPT = "experiments/t3_decomp/build_k24_form.py"
SCRIPT_PATH = Path(__file__).resolve()
BRANCH = bo.BRANCH
FROZEN_FACTOR = bo.FROZEN_FACTOR
TABLE_M, TABLE_SHA256 = bo.TABLE_M, bo.TABLE_SHA256
BASES = bo.BASES
KINDS = {"k24sc3": ("x_sc3", "sc3", ("sc3",)), "k24bulk": ("x_bulk", "bulk", ("bulk",)),
         "k24joint": ("x_joint", "joint", ("sc3", "bulk"))}  # form kind -> (K24 arm, k24_table kind, table keys)
FORMS = tuple(f"{k}-{b}" for k in KINDS for b in BASES)
READING = HERE / "results" / "k24_wider_tables" / "reading.json"
REGISTER = REPO / "EXTERNAL_DATA.md"
STUDY_SERIES = {"GSE204908": ["GSE204906"], "GSE277839": ["GSE277839"], "GSE299068": ["GSE299068"],
                "MGP-DMDD": ["E-MTAB-7322"]}  # study -> the accession its register row names
PRICED_AT = "x = 1 only"
LIST_ENTRY = "experiments/harness/TRIES.md row K24 and decision log 2026-11-21 (c)"
EXIT_REFUSED, EXIT_FORM_REFUSED = 2, 4

Refused = bs.Refused


class FormRefused(Refused):
    """The form's own refusal (exit status 4)."""


app = typer.Typer(no_args_is_help=True, add_completion=False)


def parse_form(form) -> tuple[str, str]:
    f = str(form) if form is not None else ""
    if f not in FORMS:
        raise Refused(f"--form {form!r} is not one of {list(FORMS)} (controls, twins, x_e2, barred_form and identity "
                      "are never built here)")
    kind, base = f.split("-", 1)
    return kind, base


def register_rows(studies: list[str], path: Path | None = None) -> dict:
    """{study: [register row ids]} from EXTERNAL_DATA.md's rows naming the study's accession; Refused when none."""
    path = REGISTER if path is None else Path(path)
    rows = [ln for ln in Path(path).read_text().splitlines() if re.match(r"\| X\d+ ", ln)]
    out = {}
    for s in studies:
        accs = STUDY_SERIES.get(s, [s])
        ids = sorted({ln.split(" ", 2)[1] for ln in rows if any(a in ln for a in accs)})
        if not ids:
            raise Refused(f"EXTERNAL_DATA.md holds no row naming {accs} (study {s}): not built")
        out[s] = {"accessions": accs, "rows": ids}
    return out


def check_reading(arm: str, path: Path | None = None, check_tracked: bool = True) -> tuple[dict, dict]:
    """The committed K24 reading names exactly the arm, is not VOID or incomplete, records the listing and pins and
    holds the arm undropped; returns (its record, its pre-gate)."""
    path = READING if path is None else Path(path)
    if path.is_symlink() or not path.is_file():
        raise Refused(f"no K24 reading at {path.name}: a K24 form is built only when the reading names it")
    if check_tracked and not ext.tracked(path):
        raise Refused(f"{path.name} is not committed: a K24 form is built only from a committed reading")
    raw = path.read_bytes()
    r = json.loads(raw)
    if r.get("incomplete") or r.get("void") is not False:
        raise Refused("the K24 reading is VOID or incomplete: it names no file")
    if (r.get("names") or []) != [arm]:
        raise Refused(f"the K24 reading does not name exactly {arm} (it names {r.get('names') or 'none'}): not built")
    e = r.get("ext") or {}
    if e.get("listing_sha256") != k24.K24_LISTING_SHA256 or e.get("tables") != dict(k24.K24_TABLE_SHA256):
        raise Refused("the K24 reading records other tables than the listing and pins: not built")
    pre = r.get("pregate") or {}
    if arm not in k24.active_arms(pre) or arm not in (r.get("arms") or []):
        raise Refused(f"the K24 reading's pre-gate holds {arm} dropped: not built")
    return ({"file": str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path),
             "sha256": hashlib.sha256(raw).hexdigest(), "names": r["names"], "outcome": r.get("outcome")}, pre)


def consts_from(pre: dict, kind: str, tables: dict) -> dict:
    """e_T and c from the reading's pre-gate, each recomputed from the loaded tables to a relative 1e-12."""
    out = {}
    for key in ("sc3", "bulk"):
        if key not in tables:
            continue
        want = (pre.get("e") or {}).get(key)
        got = k24.exponent(tables[key])
        if want is None or got is None or not math.isclose(float(want), got, rel_tol=k24.REL_TOL, abs_tol=0.0):
            raise Refused(f"e_{key}: the reading's {want!r} is not the one the loaded table gives ({got!r}): not "
                          "built")
        out[f"e_{key}"] = float(want)
    if kind == "joint":
        want = pre.get("c")
        got = k24.joint_c(tables["sc3"], tables["bulk"], out["e_sc3"], out["e_bulk"])
        if want is None or got is None or not math.isclose(float(want), got, rel_tol=k24.REL_TOL, abs_tol=0.0):
            raise Refused(f"c: the reading's {want!r} is not the one the loaded tables give ({got!r}): not built")
        out["c"] = float(want)
    return out


def form_factor(kind: str, x: float, phi_m: np.ndarray, tables: dict, consts: dict) -> tuple[np.ndarray, dict]:
    """phi_M ** x * k24_table(kind); FormRefused where the product cannot be made or changes fewer than MIN_CHANGED
    genes against phi_M ** x."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    base = phi_m if x == 1.0 else phi_m ** x
    try:
        scaled = k24.k24_table(kind, tables, consts)
        phi = ext.product(phi_m, x, scaled)
        figs = ext.figures(base, phi)
    except ext.FormRefused as e:
        raise FormRefused(f"{kind}: {e}") from e
    if figs["genes_changed"] < ext.MIN_CHANGED:
        raise FormRefused(f"{kind}: the product changes {figs['genes_changed']} genes against phi_M ** x "
                          f"(< {ext.MIN_CHANGED})")
    return phi, {"applied": figs, "scaled": ext.table_figures(phi_m, scaled)}


@contextmanager
def k24_loader(kind: str, x: float, pre: dict, out: dict, env: dict | None = None, register: Path | None = None):
    """bs.load_factor, for this process only: the frozen loader for phi_M (checked against TABLE_SHA256 too), then the
    form's tables through k24_wider_tables (by key, against the same genes) and their register rows, then the form's
    factor; every refusal comes before any cell is drawn."""
    original = bs.load_factor
    tkind, keys = KINDS[kind][1], KINDS[kind][2]

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        phi_m, rec = original(run_dir, TABLE_M, genes, dec)
        if rec["sha256"] != TABLE_SHA256:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[:16]}): not built")
        try:
            listing, lsha = k24.read_listing()
            loaded = {k: k24.load_k24_table(k, list(genes), listing, lsha, env=env) for k in keys}
        except k24.Refused as e:
            raise Refused(str(e)) from None
        tables = {k: v[0] for k, v in loaded.items()}
        reg = register_rows(sorted({s for v in loaded.values() for s in v[1]["studies"]}), register)
        consts = consts_from(pre, tkind, tables)
        phi, frec = form_factor(tkind, x, phi_m, tables, consts)
        out.update(table_m=rec, tables={k: v[1] for k, v in loaded.items()}, consts=consts, factor=frec, register=reg,
                   n_loads=out.get("n_loads", 0) + 1)
        return phi, rec

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


TABLE_FIELDS = ("key", "name", "file", "dir_realpath", "sha256", "listing_sha256", "pinned_sha256", "combine_file",
                "combine_sha256", "manifest_sha256", "modality", "studies", "tiers", "sign_rule", "refstats_sha256")


def build_form(decision: Path, form, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None = None,
               seed: int = bo.LIST_SEED, overwrite: bool = False, threads: int = 16, env: dict | None = None,
               reading: Path | None = None, register: Path | None = None) -> dict:
    kind, base = parse_form(form)
    slot, x, base_sha = BASES[base]
    frozen = bo.check_frozen()
    tree = bo.check_clean()
    try:
        k24.check_pins(k24.K24_LISTING_SHA256, k24.K24_TABLE_SHA256)
    except k24.Refused as e:
        raise Refused(str(e)) from None
    arm = KINDS[kind][0]
    rrec, pre = check_reading(arm, reading)
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    got: dict = {}
    with k24_loader(kind, x, pre, got, env, register):
        side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                   run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the form's factor once at w = 1")
    side["inputs"]["factor_note"] = bo.FACTOR_NOTE
    reg = got["register"]
    side["outside_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH),
        "k24_wider_tables_sha256": bs.sha256(HERE / "k24_wider_tables.py"),
        "outside_forms_sha256": bs.sha256(HERE / "outside_forms.py"), "form": f"{kind}-{base}", "arm": arm,
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha},
        "table_m": {"file": got["table_m"]["file"], "sha256": got["table_m"]["sha256"], "pinned_sha256": TABLE_SHA256},
        "tables": {k: {**{f: r.get(f) for f in TABLE_FIELDS},
                       "register": {s: reg[s] for s in r["studies"]}} for k, r in got["tables"].items()},
        "e": {k: v for k, v in got["consts"].items() if k.startswith("e_")}, "c": got["consts"].get("c"),
        "s_star": k24.S_STAR, "factor": got["factor"], "frozen_files_sha256": frozen, "tree_at_start": tree,
        "k24_reading": rrec, "priced_at": PRICED_AT, "list_entry": LIST_ENTRY, "depth_step": None,
    }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


@app.callback()
def main() -> None:
    """The barred file with its Mab21l2 factor times one of row K24's scaled outside-knockout tables."""


@app.command()
def build(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    form: str = typer.Option(..., help=f"one of {', '.join(FORMS)}; nothing else"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the frozen builder"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    out: Path = typer.Option(..., help="output .h5ad (sidecar: same name with .json)"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(bo.LIST_SEED, help="1000001, the list's seed; any other is refused"),
    n_cells: str = typer.Option(bo.LIST_N_CELLS, help="'max', the list's count; any other is refused"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    overwrite: bool = typer.Option(False, help="replace an existing output"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """Build one fixed form."""
    try:
        if seed != bo.LIST_SEED or n_cells != bo.LIST_N_CELLS:
            raise Refused(bo.LIST_INPUTS_REFUSAL)
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
    ob = s["outside_build"]
    bs.say(f"form {ob['form']}: base {ob['base']['name']} (x = {ob['base']['x']}), genes changed against phi_M ** x "
           f"{ob['factor']['applied']['genes_changed']}")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


if __name__ == "__main__":
    app()
