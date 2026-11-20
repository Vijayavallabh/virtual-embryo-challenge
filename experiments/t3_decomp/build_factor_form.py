#!/usr/bin/env python
"""The barred branch's T3 file with its Mab21l2 factor in one of the regularised forms of try K21 (ledger row K21 and
the decision-log entry that fixes its names before the run). A file is built only in a form the K21 reading names, on
the base B that decision log 2026-11-20 (q) fixes; nothing here chooses a form or a base.

A thin builder, in the manner of build_two_knockout.py (which it does not change). It imports the frozen builder
(build_submission.py) and script (t3_decompose.py) unchanged and builds a barred-branch file exactly as
`build_submission.py build DECISIONS --branch barred --submission s1|s2` does, with one change, for this process only:

  top100, cap90   the factor the frozen loader returns for `transfer_mab21l2_all` (phi_M, sha256-checked by the frozen
                  loader against the decision file and here against the pin) is replaced by T(phi_M) ** x, T the
                  form's transform in factor_forms.py (the module try K21's arms apply), its constants fixed from phi_M
                  itself; nothing is rescaled or clipped afterwards; t3_decompose.build receives it at w = 1
  wtprod          the frozen construct is called with no factor (the same rows: build draws every row before its factor
                  step), then the frozen apply_wt_revert moves each cell by lineage group at exponent 0.5 with
                  per'_g = (phi_M ** x)**2 * per_g and the fallback (phi_M ** x)**2 * phi_pooled, so each
                  group's applied factor is phi_M ** x * per_g ** 0.5; per_g and phi_pooled from the frozen
                  wt_revert_factors on the six other released WT stage files, on the stream
                  "submission:barred:<slot>:wt_revert"
  identity        phi_M ** x: the base file byte for byte (BASES' sha256), the build launcher's gate

Forms (--form; nothing else is accepted, no number is taken from the command line): identity-B, top100-B, cap90-B and
wtprod-B for B in s1, s1x1.5, s1x2.0, s2. The controls (uni_top, uni_cap, wtprod_perm) and the dropped forms (soft50s,
rand100) are refused by name. top100 is refused when 100 or more genes sit at the clip or |ln phi| ties across the cut.

Refusals (exit status 2), before anything loads: another form; a seed other than 1000001 or an --n-cells other than
'max' on the command line; frozen files whose sha256 differs from the pins (the same bytes build_two_knockout.py and
build_barred_exponent.py pin); a tree that is not clean (vec.provenance's rule); the frozen builder's decision-file
check with no allow flag; when the factor loads, a table whose sha256 differs from the decision file's record or from
the pin, or whose genes differ from the reference's; a base with a depth, stage or coverage step (wtprod).

The sidecar is the frozen builder's plus one block `form_build` (this script and its sha256, the form, the base, x, the
transform's record, the factor's figures, the table and its hash, the WT series hashes for wtprod, the pinned frozen
hashes, the ledger row). It adds no cell number. `inputs.factor` stays the frozen loader's record of phi_M; the factor
applied is the one `form_build` names, as `inputs.factor_note` says.

  uv run --frozen --no-sync python -I -B experiments/t3_decomp/build_factor_form.py build DECISIONS.json \\
      --form top100-s1 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> --board T3:bcatenin \\
      --seed 1000001 --out runs/<folder>/T3_top100_s1.h5ad
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
import build_submission as bs  # noqa: E402
import factor_forms as ff  # noqa: E402

from vec import provenance  # noqa: E402

SCRIPT = "experiments/t3_decomp/build_factor_form.py"
SCRIPT_PATH = Path(__file__).resolve()
BRANCH = "barred"
FROZEN_FACTOR = {"source": "mab21l2", "exponent": 1.0}
TABLE_M = "transfer_mab21l2_all"
TABLE_SHA256 = "64b1e0b8572c7193ed1e1f0b7339a26081e601f1e706e6b08193e1718dc5c7e8"  # gata4_decisions.json's record
# the frozen files' bytes, as build_barred_exponent.py and build_two_knockout.py pin them (a test checks the equality)
FROZEN_SHA256 = {
    "t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
    "rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13",
    "build_submission.py": "1194a63993a27dde4dedfdb9180de035c5f8f2a57aaa7d97b38e7a15ea9e0c74",
}
# base -> (frozen slot, x, the reference sha256 of the base file built on the DGX), as build_two_knockout.BASES
BASES = {
    "s1": ("s1", 1.0, "3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2"),
    "s1x1.5": ("s1", 1.5, "c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420"),
    "s1x2.0": ("s1", 2.0, "7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661"),
    "s2": ("s2", 1.0, "1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f"),
}
KINDS = ("identity", "top100", "cap90", "wtprod")
FORMS = tuple(f"{k}-{b}" for k in KINDS for b in BASES)
REFUSED_KINDS = {"uni_top": "a control of try K21, never built", "uni_cap": "a control of try K21, never built",
                 "wtprod_perm": "a control of try K21, never built",
                 "soft50s": "dropped before try K21 (its size-fixing constant sends the clip genes past 4)",
                 "rand100": "dropped before try K21 (not size-matched)", "barred_form": "use identity-<base>"}
LIST_ENTRY = "experiments/harness/TRIES.md row K21 and the decision-log entry that fixes its names"
LIST_SEED, LIST_N_CELLS = 1_000_001, "max"
LIST_INPUTS_REFUSAL = "the list fixes seed 1000001 and max_cells"
FACTOR_NOTE = "the frozen loader's record of phi_M only; the factor applied is form_build.factor"
WT_REVERT_STREAM = bs.WT_REVERT_STREAM  # "wt_revert": the frozen builder's own stream name

Refused = bs.Refused
app = typer.Typer(no_args_is_help=True, add_completion=False)


def parse_form(form) -> tuple[str, str]:
    """--form: exactly one of FORMS; returns (kind, base), else Refused (controls and dropped forms by name)."""
    f = str(form) if form is not None else ""
    kind = f.split("-", 1)[0]
    if kind in REFUSED_KINDS and f not in FORMS:
        raise Refused(f"--form {form!r} is refused: {kind} is {REFUSED_KINDS[kind]}")
    if f not in FORMS:
        raise Refused(f"--form {form!r} is not one of {list(FORMS)}")
    kind, base = f.split("-", 1)
    return kind, base


def frozen_paths() -> dict[str, Path]:
    return {"t3_decompose.py": bs.T3D_PATH, "rules.yaml": bs.RULES_PATH, "build_submission.py": bs.BUILDER_PATH}


def check_frozen(paths: dict[str, Path] | None = None) -> dict[str, str]:
    paths = frozen_paths() if paths is None else paths
    got = {k: bs.sha256(p) for k, p in paths.items()}
    bad = [k for k in FROZEN_SHA256 if got.get(k) != FROZEN_SHA256[k]]
    if bad:
        raise Refused(f"frozen file(s) {bad} differ from the pinned bytes: not built")
    return got


def check_clean() -> dict:
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    if commit == "unknown" or unclean is None:
        raise Refused("git cannot tell the tree's state: not built (the stamp would not say clean)")
    if unclean:
        raise Refused(f"the tree is not clean ({len(unclean)} path(s), e.g. {unclean[:3]}): commit first (the build "
                      "stamp must say clean)")
    return {"commit": commit, "clean": True}


def form_factor(kind: str, x: float, phi_m: np.ndarray) -> tuple[np.ndarray, dict]:
    """The vector applied for a vector form: T(phi_M) ** x (factor_forms.vector_form, then at_exponent); for wtprod,
    phi_M ** x (the construct step applies the product). Raises Refused where the form cannot be made."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    try:
        name = "identity" if kind == "wtprod" else kind
        t, rec = ff.vector_form(name, phi_m)
    except ff.FormRefused as e:
        raise Refused(f"{kind}: {e}") from e
    phi = ff.at_exponent(t, x)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the form's factor is not finite and positive")
    return phi, {"transform": rec, "transformed_at_x1": ff.describe(phi_m, t), "applied": ff.describe(phi_m, phi)}


@contextmanager
def form_loader(kind: str, x: float, out: dict):
    """bs.load_factor, for this process only: the frozen loader for phi_M, checked against TABLE_SHA256 too, then the
    form's factor."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        phi_m, rec = original(run_dir, TABLE_M, genes, dec)
        if rec["sha256"] != TABLE_SHA256:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[:16]}): not built")
        phi, frec = form_factor(kind, x, phi_m)
        out.update(table=rec, factor=frec, phi_m=phi_m, n_loads=out.get("n_loads", 0) + 1)
        return phi, rec

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


@contextmanager
def wtprod_construct(wt: Path, series: dict, seed: int, slot: str, out: dict):
    """bs.construct, for this process only: the frozen construct with no factor (the same rows), then the frozen
    apply_wt_revert with per'_g = phi**2 * per_g and the fallback phi**2 * phi_pooled (phi = phi_M ** x, as the loader
    returned it), on the frozen builder's stream name for the step."""
    original = bs.construct

    def construct(ref, n, rng, f, phi=None, w=1.0, ratio=None, alphas=None):
        if phi is None or w != 1.0:
            raise Refused("wtprod expects the barred factor at w = 1")
        if ratio is not None or alphas is not None:
            raise Refused("wtprod is for a base with no depth, stage or coverage step")
        built = original(ref, n, rng, f, phi=None, w=1.0, ratio=None, alphas=None)
        pooled, per, group, rrec = bs.wt_revert_factors(Path(wt), series, ref.genes, ref.pool.C)
        per2, fb2 = ff.wtprod_factors(phi, per, pooled)
        stream = f"submission:{BRANCH}:{slot}:{WT_REVERT_STREAM}"
        moved, how = bs.apply_wt_revert(ref, built, fb2, per2, group, bs.t3d.stable_rng(seed, stream))
        applied = {g: v ** ff.WT_EXPONENT for g, v in {**per2, "pooled_fallback": fb2}.items()}
        out.update(wt_revert={k: v for k, v in rrec.items()}, wt_stream=stream, by_group=how["factor_by_group"],
                   applied={"min": float(min(v.min() for v in applied.values())),
                            "max": float(max(v.max() for v in applied.values())),
                            "sd_ln_by_group": {g: ff.sd(np.log(v)) for g, v in applied.items()}},
                   n_constructs=out.get("n_constructs", 0) + 1)
        return moved

    bs.construct = construct
    try:
        yield
    finally:
        bs.construct = original


def build_form(decision: Path, form, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None = None,
               seed: int = 1_000_001, overwrite: bool = False, threads: int = 16,
               wt_series: dict[str, Path] | None = None) -> dict:
    """The named form through bs.build_submission; returns the sidecar (the frozen builder's plus `form_build`)."""
    kind, base = parse_form(form)
    slot, x, base_sha = BASES[base]
    frozen = check_frozen()
    tree = check_clean()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    series = None
    if kind == "wtprod":
        if any(spec.get(k) for k in bs.EXTRAS):
            raise Refused(f"wtprod: the barred {slot} carries a step {[k for k in bs.EXTRAS if spec.get(k)]}")
        series = bs.wt_revert_series(Path(wt), wt_series)
    got: dict = {}
    wgot: dict = {}
    with form_loader(kind, x, got):
        if kind == "wtprod":
            with wtprod_construct(Path(wt), series, seed, slot, wgot):
                side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                           run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        else:
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR or (kind == "wtprod" and wgot.get("n_constructs") != 1):
        raise RuntimeError("the build did not apply the form's factor once at w = 1")
    side["inputs"]["factor_note"] = FACTOR_NOTE
    factor = got["factor"]
    if kind == "wtprod":
        factor = factor | {"wtprod": {"applied": "phi_M ** x * per_g ** 0.5 per lineage group (fallback phi_M ** x * "
                                                 "phi_pooled ** 0.5)", "stream": wgot["wt_stream"],
                                      "factor_by_group": wgot["by_group"], "applied_range": wgot["applied"],
                                      "wt_revert": wgot["wt_revert"]}}
        side["inputs"]["rng_streams"] = [*side["inputs"]["rng_streams"], wgot["wt_stream"]]
    side["form_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH),
        "factor_forms": "experiments/t3_decomp/factor_forms.py",
        "factor_forms_sha256": bs.sha256(HERE / "factor_forms.py"), "form": f"{kind}-{base}",
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha},
        "constants": {"top_k": ff.TOP_K, "cap_q": ff.CAP_Q, "clip": list(ff.CLIP), "clip_tol": ff.CLIP_TOL,
                      "wt_exponent": ff.WT_EXPONENT},
        "table": {"file": got["table"]["file"], "sha256": got["table"]["sha256"], "pinned_sha256": TABLE_SHA256},
        "factor": factor, "frozen_files_sha256": frozen, "tree_at_start": tree, "list_entry": LIST_ENTRY,
        "identity": f"the {base} base file byte for byte ({base_sha[:12]})" if kind == "identity" else None,
        "depth_step": None,
    }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


@app.callback()
def main() -> None:
    """The barred file with its Mab21l2 factor in one of try K21's fixed forms."""


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
    wt_series: list[str] = typer.Option(None, help="wtprod only: NAME=PATH for the six other WT stages, as the frozen "
                                        "--wt-series (default: data/raw/T2/<NAME>.h5ad)"),
):
    """Build one fixed form."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        kind, _ = parse_form(form)
        if wt_series and kind != "wtprod":
            raise Refused("--wt-series is for wtprod only")
        b = bs.board_spec(board, boards_dir)
        s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                       overwrite=overwrite, threads=threads,
                       wt_series=bs.parse_pairs(wt_series, "--wt-series") if wt_series else None)
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    chk = s["output"]["check"]
    fb = s["form_build"]
    bs.say(f"form {fb['form']}: base {fb['base']['name']} (x = {fb['base']['x']}), genes changed against phi_M "
           f"{fb['factor']['applied']['genes_changed']}")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


if __name__ == "__main__":
    app()
