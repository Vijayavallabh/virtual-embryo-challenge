#!/usr/bin/env python
"""The barred branch's T3 file with its Mab21l2 factor multiplied by the released Gata4 knockout's factor: the forms of
the two-knockout list (ledger row R17; the decision-log entry that fixes the list before any of its scores).

A thin builder, in the manner of build_barred_exponent.py. It imports the frozen builder (build_submission.py) and
script (t3_decompose.py) unchanged and builds a barred-branch file exactly as `build_submission.py build DECISIONS
--branch barred --submission s1|s2` does, with one change: the factor the frozen loader returns for the base's table is
replaced, for this process only, by a fixed combination of two committed factor tables, both loaded and sha256-checked
by the frozen loader (`build_submission.load_factor`) and checked again against the sha256 pinned below:

  phi_M  <run>/factors/transfer_mab21l2_all.csv   the released Mab21l2 knockout against WT E8.75 (the barred factor)
  phi_G  <run>/factors/gata4_pooled.csv            the released Gata4 knockout (both validation replicates) against
                                                   WT E8.75, strict non-Mesp1 labels (the allowed branch's factor)

Each table is used as loaded: t3_decompose.gene_factor already rescaled it (den-weighted total proportion unchanged)
and clipped it to [0.25, 4] when the run wrote it. After the combination nothing is rescaled or clipped again, as the
exponent builder (row R14) and the frozen build treat every exponent. So the combined factor is not rescaled: its
common shift is recorded in the sidecar (gene mean of ln phi, cos^2 of ln phi with the all-ones vector), not treated.

The forms (--form; nothing else is accepted, and no number is taken from the command line):

  identity-s1        phi_M                           = the frozen barred s1 byte for byte (3712c3f9988d...)
  identity-s1x1.5    phi_M ** 1.5                    = build_barred_exponent.py --exponent 1.5 (c060fcaa6c7f...)
  identity-s1x2.0    phi_M ** 2.0                    = build_barred_exponent.py --exponent 2.0 (7393810fbe3b...)
  identity-s2        phi_M                           = the frozen barred s2 (1b69521e3c28...)
  c1-s1, c1-s2       phi_M * phi_G                   (C1: base x = 1 only)
  c2-s1, c2-s1x1.5,  phi_M ** (x c) * phi_G ** c     (C2: held at the base's size)
  c2-s1x2.0, c2-s2   c = SD(x ln phi_M) / SD(x ln phi_M + ln phi_G), population SD over the panel's genes,
                     computed here from the two tables as loaded, never from a score

C1 on the exponent-1.5 and exponent-2.0 bases (c1-s1x1.5, c1-s1x2.0) is refused: the list fixes C1 on x = 1 only.

Where each quantity enters. The base fixes the slot (s1: the central composition at m = D7m, cells drawn on stream
"submission:barred:s1"; s2: the mild composition, stream "submission:barred:s2"), the "...:assemble" stream, the
seed and the board; all are the frozen builder's own, so the same wild-type cells are drawn in every form of a base.
The spec keeps the frozen factor {mab21l2, exponent 1.0}, so t3_decompose.build receives w = 1 and applies the
combined phi as given (`apply_gene_factor`, after the cells are drawn, before lognorm). The base's exponent x enters
as `phi_M ** x` here, the same numpy operation the frozen build applies for w = x, so each identity form gives its
base's bytes. No depth, stage or coverage step: the barred spec has none.

Refusals (exit status 2), before anything loads: another form (taken as typed, no whitespace stripped); on the command
line a seed other than 1000001 or an --n-cells other than 'max' (the list fixes both); t3_decompose.py, rules.yaml
or build_submission.py whose sha256 differs from the bytes pinned below (the same bytes build_barred_exponent.py
pins); a tree that is not clean by vec.provenance's rule; the frozen builder's decision-file check, with no allow
flag; then, when the factor loads, a table whose sha256 differs from the decision file's record (the frozen loader) or
from the value pinned here, a table whose genes differ from the reference's (order included; the frozen loader), or a
combination that is not finite and positive.

The sidecar is the frozen builder's plus one block `combined_build` (this script and its sha256, the form, the base,
x, the Gata4 weight, c, both tables with their paths and hashes, the SDs and correlation, the common-shift figures, the
pinned frozen hashes, the list's entry). It adds no cell number. `inputs.factor` stays the frozen loader's record of
phi_M; the factor applied is the one `combined_build` names, as `inputs.factor_note` says in every form.

  uv run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build DECISIONS.json \\
      --form c1-s1 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> --board T3:bcatenin \\
      --seed 1000001 --out runs/<folder>/T3_c1_s1.h5ad

--gata4-exponent 2.0 --captain-ruling R22 (row R22's second T3 search line, the Gata4-table weight w on C1, as
amended by row R30; off by default: without both, every form, byte and behaviour above is unchanged and nothing below
is read or computed): C1 on the barred s2 with its Gata4 weight raised from 1 to 2, the one file T3_c1_s2_g2.h5ad,

  c1-s2 --gata4-exponent 2.0    phi_M ** 1 * phi_G ** 2.0     nothing rescaled or clipped after the product

applied where C1's factor is applied (bs.load_factor's replacement, at w = 1, once, after the cells are drawn).
Refused (exit 2, nothing built): any other value; any other form; --gata4-exponent without --captain-ruling R22;
--captain-ruling R22 without --gata4-exponent 2.0; any other ruling; and the pair while experiments/harness/TRIES.md
does not hold exactly one row R22 whose first cell passes build_outside_form.captain_order_names for the option
"`--gata4-exponent 2.0`". Before the build, from released WT E8.75 (den as build_outside_form.wt_den makes it) and
the two pinned tables only, the line's figures are computed and recorded: the den-weighted totals of phi_M, of
phi_M * phi_G (C1's) and of phi_M * phi_G ** 2, and the two relative shifts against phi_M's. By row R30 the line's old
5% limit gates nothing: no value of a shift stops the build (the old value is recorded as lifted). The factor the build
applies must be the very array those figures were computed from (its sha256, over WT E8.75's genes); if it is not,
exit 5 with nothing written. The sidecar's combined_build then also holds `gata4_line` (the ruling row, the line, the
weight, row R30's lifting, the figures and the sha256 of the factor applied).
"""

from __future__ import annotations

import hashlib
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

from vec import provenance  # noqa: E402

SCRIPT = "experiments/t3_decomp/build_two_knockout.py"
SCRIPT_PATH = Path(__file__).resolve()
BRANCH = "barred"
FROZEN_FACTOR = {"source": "mab21l2", "exponent": 1.0}
TABLE_M, TABLE_G = "transfer_mab21l2_all", "gata4_pooled"
# the two committed factor tables of the T3 run gata4_20261028T122425Z, as the committed decision file records them
# (experiments/t3_decomp/results/gata4/gata4_decisions.json, submissions.parameters.factor_files)
TABLE_SHA256 = {
    TABLE_M: "64b1e0b8572c7193ed1e1f0b7339a26081e601f1e706e6b08193e1718dc5c7e8",
    TABLE_G: "36f2d65a36d2d4fd5307cddfadb56e9724c362c982795f30428eda752b4734ef",
}
# the frozen files' bytes on main at ca2ebc5 (unchanged since the third rehearsal, dc7b6ae), as build_barred_exponent
FROZEN_SHA256 = {
    "t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
    "rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13",
    "build_submission.py": "1194a63993a27dde4dedfdb9180de035c5f8f2a57aaa7d97b38e7a15ea9e0c74",
}
# base -> (frozen slot, x, the reference sha256 of the base file built on the DGX)
BASES = {
    "s1": ("s1", 1.0, "3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2"),
    "s1x1.5": ("s1", 1.5, "c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420"),
    "s1x2.0": ("s1", 2.0, "7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661"),
    "s2": ("s2", 1.0, "1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f"),
}
# form kind -> the Gata4 table's weight before c (c is computed for c2 only, 1 otherwise)
KINDS = {"identity": 0.0, "c1": 1.0, "c2": 1.0}
FORMS = tuple([f"identity-{b}" for b in BASES] + ["c1-s1", "c1-s2"] + [f"c2-{b}" for b in BASES])
REFUSED_FORMS = {"c1-s1x1.5": "C1 is fixed on the x = 1 bases only (barred s1, barred s2)",
                 "c1-s1x2.0": "C1 is fixed on the x = 1 bases only (barred s1, barred s2)"}
WHERE = ("the frozen loader's phi_M and phi_G (each already rescaled and clipped to [0.25, 4] by gene_factor) are "
         "combined gene by gene here; t3_decompose.build receives the combination at w = 1 and applies it through "
         "apply_gene_factor after the cells are drawn; nothing is rescaled or clipped after the combination")
LIST_ENTRY = "docs/PLAN.md decision log 2026-11-20 (q); experiments/harness/TRIES.md row R17"
LIST_SEED, LIST_N_CELLS = 1_000_001, "max"
LIST_INPUTS_REFUSAL = "the list fixes seed 1000001 and max_cells"
FACTOR_NOTE = "the frozen loader's record of phi_M only; the factor applied is combined_build.combination"
DEN_NOTE = ("not recorded: the factor files carry gene and factor only (t3_decompose.submission_parameters), so the "
            "den-weighted mean cannot be computed from the tables alone")

# Row R22's second T3 search line (the Gata4-table weight w on C1; decision log 2026-11-23 (ff), (hh)), as row R30
# (the captain's order of 2026-11-24, decision log 2026-11-24 (d)) amends its guard: the 5% limit stops nothing.
TRIES_FILE = HERE.parents[1] / "experiments" / "harness" / "TRIES.md"
G2_ROW = "R22"
G2_EXPONENTS = (2.0,)  # --gata4-exponent: the single value 2.0
G2_FORM = "c1-s2"
G2_OPTION = "`--gata4-exponent 2.0`"
G2_FILE = "T3_c1_s2_g2.h5ad"
G2_LINE = ("T3 search line 2 (R22): the Gata4-table weight w in C1's factor phi_M ** 1 * phi_G ** w on the barred s2; "
           "B (w = 0), C1 (w = 1), this file T3_c1_s2_g2.h5ad (w = 2.0); the line ends at 2.0")
G2_R30 = ("row R30 (the captain's order of 2026-11-24; decision log 2026-11-24 (d)): the limit of 5% on the "
          "den-weighted total gates nothing for this file; the figures are computed and recorded and the file is "
          "built whatever they are")
G2_OLD_LIMIT = 0.05  # the line's old guard limit; recorded only, lifted by row R30, gates nothing
G2_LIST_ENTRY = ("experiments/harness/TRIES.md rows R22 and R30; decision log 2026-11-23 (ff), (hh) and "
                 "2026-11-24 (d)")
G2_DEN_RULE = ("WT E8.75's mean within-cell proportions (build_outside_form.wt_den: load_stage, duplicates dropped); "
               "a den-weighted total is den . phi / sum den")
EXIT_REFUSED, EXIT_FACTOR_MISMATCH = 2, 5

Refused = bs.Refused


class FactorMismatch(Refused):
    """The factor the build applies is not the array the line's figures were computed from (exit status 5)."""


app = typer.Typer(no_args_is_help=True, add_completion=False)


def parse_form(form) -> tuple[str, str]:
    """--form: exactly one of FORMS; returns (kind, base), else Refused."""
    f = str(form) if form is not None else ""
    if f in REFUSED_FORMS:
        raise Refused(f"--form {form!r} is refused: {REFUSED_FORMS[f]}")
    if f not in FORMS:
        raise Refused(f"--form {form!r} is not one of {list(FORMS)}")
    kind, base = f.split("-", 1)
    return kind, base


def frozen_paths() -> dict[str, Path]:
    """The three frozen files as the frozen builder itself names them."""
    return {"t3_decompose.py": bs.T3D_PATH, "rules.yaml": bs.RULES_PATH, "build_submission.py": bs.BUILDER_PATH}


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


def _sd(v: np.ndarray) -> float:
    return float(np.std(v, ddof=0))


def _shift(phi: np.ndarray) -> dict:
    """The common shift of a factor, from the factor alone: gene mean of ln phi and cos^2(ln phi, 1)."""
    lp = np.log(phi)
    ss = float(lp @ lp)
    return {"gene_mean_ln": float(lp.mean()), "cos2_ln_with_ones": float(lp.sum() ** 2 / (lp.size * ss)) if ss > 0
            else 0.0, "sd_ln": _sd(lp)}


def combine(kind: str, x: float, phi_m: np.ndarray, phi_g: np.ndarray) -> tuple[np.ndarray, dict]:
    """The form's factor from the two tables as loaded, gene by gene; returns (phi, record).

    identity: phi_M ** x (phi_M itself at x = 1); c1: phi_M ** x * phi_G; c2: phi_M ** (x c) * phi_G ** c with
    c = SD(x ln phi_M) / SD(x ln phi_M + ln phi_G). `**` is the operation the frozen build applies for w = x."""
    if kind not in KINDS:
        raise Refused(f"unknown form kind {kind!r}")
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    if phi_m.shape != phi_g.shape or phi_m.ndim != 1:
        raise Refused("the two factor tables do not cover the same genes")
    for nm, p in ((TABLE_M, phi_m), (TABLE_G, phi_g)):
        if not np.all(np.isfinite(p) & (p > 0)):
            raise Refused(f"{nm}: factor values must be finite and positive")
    base = phi_m if x == 1.0 else phi_m ** x
    g = KINDS[kind]
    lm, lg = np.log(phi_m), np.log(phi_g)
    sd_base, sd_sum = _sd(x * lm), _sd(x * lm + lg)
    c = 1.0
    if kind == "identity":
        phi = base
    elif kind == "c1":
        phi = base * phi_g
    else:
        if not sd_sum > 0:
            raise Refused("SD(x ln phi_M + ln phi_G) is 0: c is undefined")
        c = sd_base / sd_sum
        phi = (phi_m ** (x * c)) * (phi_g ** c)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    sd_m, sd_g = _sd(lm), _sd(lg)
    rec = {"kind": kind, "x": float(x), "gata4_weight": g * c, "mab21l2_exponent": float(x) * c, "c": float(c),
           "c_rule": "SD(x ln phi_M) / SD(x ln phi_M + ln phi_G), population SD over the genes" if kind == "c2"
           else "c = 1",
           "sd_ln_phi_m": sd_m, "sd_ln_phi_g": sd_g, "sd_ln_base": sd_base, "sd_ln_base_plus_g": sd_sum,
           "corr_ln_m_ln_g": float(np.corrcoef(lm, lg)[0, 1]) if sd_m > 0 and sd_g > 0 else None,
           "base_factor": _shift(base), "factor": _shift(phi), "den_weighted_mean": DEN_NOTE,
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


@contextmanager
def combined_loader(kind: str, x: float, out: dict):
    """bs.load_factor, for this process only: the frozen loader for phi_M (the barred spec's table) and for phi_G,
    each checked against TABLE_SHA256 too, then `combine`. `out` receives both tables' records and the combination."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = combine(kind, x, phi_m, phi_g)
        out.update(tables={TABLE_M: rec_m, TABLE_G: rec_g}, combination=crec, n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


# ------------------------------------------------------------------ row R22's second T3 line (off by default)


def _outside_form():
    """build_outside_form, imported only on the R22 line (its captain_order_names and wt_den are reused, not copied)."""
    import build_outside_form

    return build_outside_form


def check_gata4(form, gata4_exponent=None, captain_ruling=None, tries: Path | None = None) -> dict | None:
    """--gata4-exponent and --captain-ruling: None when neither is given (nothing is read); else the line's record,
    or Refused: any other ruling, the ruling without the option, any value but 2.0, any form but c1-s2, the option
    without the ruling, or a ledger whose row R22 does not pass build_outside_form.captain_order_names."""
    if gata4_exponent is None and captain_ruling is None:
        return None
    if captain_ruling is not None and captain_ruling != G2_ROW:
        raise Refused(f"--captain-ruling {captain_ruling!r} is refused: only row {G2_ROW} orders a file of this "
                      "builder")
    if gata4_exponent is None:
        raise Refused(f"--captain-ruling {G2_ROW} is taken only with --gata4-exponent 2.0 on --form {G2_FORM}")
    try:
        w = float(gata4_exponent)
    except (TypeError, ValueError):
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0") from None
    if isinstance(gata4_exponent, bool) or w not in G2_EXPONENTS:
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0")
    if form != G2_FORM:
        raise Refused(f"--gata4-exponent 2.0 is taken on --form {G2_FORM} only, not {form!r}")
    if captain_ruling is None:
        raise Refused(f"--gata4-exponent 2.0 is taken only with --captain-ruling {G2_ROW}")
    p = TRIES_FILE if tries is None else Path(tries)
    if not _outside_form().captain_order_names(G2_ROW, G2_OPTION, p):
        raise Refused(f"--gata4-exponent 2.0 --captain-ruling {G2_ROW} is refused: experiments/harness/TRIES.md has "
                      f"no single row {G2_ROW} whose first cell is the captain's order (\"ruling, no run\") naming "
                      f"{G2_OPTION} and `--captain-ruling {G2_ROW}`")
    return {"row": G2_ROW, "option": G2_OPTION, "gata4_exponent": w, "ledger": "experiments/harness/TRIES.md",
            "ledger_sha256": bs.sha256(p)}


def factor_sha256(phi) -> str:
    """The sha256 of a factor's float64 bytes, as build_outside_form's guard records it."""
    return hashlib.sha256(np.ascontiguousarray(phi, dtype=np.float64).tobytes()).hexdigest()


def c1_gata4_weighted(phi_m, phi_g, w: float) -> tuple[np.ndarray, dict]:
    """C1 on an x = 1 base with the Gata4 weight w: phi_M ** 1 * phi_G ** w (phi_M ** 1 is phi_M itself, as combine's
    base at x = 1), nothing rescaled or clipped after the product; returns (phi, record)."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    if phi_m.shape != phi_g.shape or phi_m.ndim != 1:
        raise Refused("the two factor tables do not cover the same genes")
    for nm, p in ((TABLE_M, phi_m), (TABLE_G, phi_g)):
        if not np.all(np.isfinite(p) & (p > 0)):
            raise Refused(f"{nm}: factor values must be finite and positive")
    w = float(w)
    phi = phi_m * (phi_g ** w)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    lm, lg = np.log(phi_m), np.log(phi_g)
    sd_m, sd_g = _sd(lm), _sd(lg)
    rec = {"kind": "c1", "x": 1.0, "gata4_weight": w, "mab21l2_exponent": 1.0, "c": 1.0,
           "c_rule": "c = 1; the Gata4 weight is the line's w", "sd_ln_phi_m": sd_m, "sd_ln_phi_g": sd_g,
           "sd_ln_base": _sd(lm), "sd_ln_base_plus_wg": _sd(lm + w * lg),
           "corr_ln_m_ln_g": float(np.corrcoef(lm, lg)[0, 1]) if sd_m > 0 and sd_g > 0 else None,
           "base_factor": _shift(phi_m), "factor": _shift(phi),
           "den_weighted_mean": "recorded in combined_build.gata4_line.figures",
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


def gata4_figures(phi_m, phi_g, w: float, den) -> dict:
    """The line's figures, a pure function of the two tables and den: the den-weighted totals of phi_M, of
    phi_M * phi_G (C1's) and of phi_M * phi_G ** w, and the relative shifts of the last two against phi_M's. They
    gate nothing (row R30); the old limit is recorded as lifted. Also the sha256 of the factor they describe."""
    phi, _ = c1_gata4_weighted(phi_m, phi_g, w)
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    if den.shape != phi_m.shape or not np.all(np.isfinite(den) & (den >= 0)) or not den.sum() > 0:
        raise Refused("den does not cover the tables' genes with finite, non-negative weights and a positive sum")
    s = float(den.sum())
    t_m, t_c1, t_w = (float(den @ f / s) for f in (phi_m, phi_m * phi_g, phi))
    return {"gata4_exponent": float(w), "den_total_phi_m": t_m, "den_total_c1": t_c1, "den_total_factor": t_w,
            "den_shift_c1": t_c1 / t_m - 1.0, "den_shift_factor": t_w / t_m - 1.0,
            "old_limit": {"max_abs_den_shift": G2_OLD_LIMIT, "status": "lifted by row R30: gates nothing"},
            "gates": None, "factor_sha256": factor_sha256(phi)}


def run_gata4_figures(dec: dict, wt: Path, run_dir: Path | None, w: float) -> dict:
    """The figures from released WT E8.75 and the two pinned tables only, before any build: den and genes by
    build_outside_form.wt_den, both tables through the frozen loader in that gene order (checked against
    TABLE_SHA256). Never refuses on a figure."""
    genes, den = _outside_form().wt_den(Path(wt))
    phi_m, rec_m = bs.load_factor(run_dir, TABLE_M, genes, dec)
    phi_g, rec_g = bs.load_factor(run_dir, TABLE_G, genes, dec)
    for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
        if rec["sha256"] != TABLE_SHA256[nm]:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[nm][:16]}): not built")
    den = np.asarray(den, dtype=np.float64)
    g = gata4_figures(phi_m, phi_g, w, den)
    g |= {"wt": str(wt), "wt_sha256": bs.sha256(Path(wt)), "den_rule": G2_DEN_RULE, "den": den.tolist(),
          "den_sha256": factor_sha256(den), "n_genes": len(genes),
          "genes_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
          "tables_sha256": {TABLE_M: rec_m["sha256"], TABLE_G: rec_g["sha256"]}, "genes": list(genes)}
    return g


@contextmanager
def gata4_loader(w: float, out: dict, figures: dict):
    """bs.load_factor, for this process only, on the R22 line: the frozen loader for phi_M and phi_G (checked against
    TABLE_SHA256 too), then c1_gata4_weighted at w, the hook C1's factor goes through. FactorMismatch (exit 5) unless
    the build's genes are the figures' and the factor's sha256 is the one the figures were computed from."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        if list(genes) != figures["genes"]:
            raise FactorMismatch("the build's genes are not those the line's figures were computed over: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = c1_gata4_weighted(phi_m, phi_g, w)
        sha = factor_sha256(phi)
        if sha != figures["factor_sha256"]:
            raise FactorMismatch(f"the factor the build would apply ({sha[:16]}) is not the one the line's figures "
                                 f"were computed from ({figures['factor_sha256'][:16]}): not built")
        out.update(tables={TABLE_M: rec_m, TABLE_G: rec_g}, combination=crec, factor_sha256=sha,
                   n_loads=out.get("n_loads", 0) + 1)
        return phi, rec_m

    bs.load_factor = load
    try:
        yield
    finally:
        bs.load_factor = original


def build_form(decision: Path, form, wt: Path, board, n_cells: int, out: Path, run_dir: Path | None = None,
               seed: int = 1_000_001, overwrite: bool = False, threads: int = 16, gata4_exponent=None,
               captain_ruling=None, tries: Path | None = None) -> dict:
    """The named form through bs.build_submission; returns the sidecar (written next to the file, the frozen
    builder's plus `combined_build`). gata4_exponent 2.0 with captain_ruling 'R22' (default None, today's paths
    unchanged) builds row R22's second T3 line's file, c1-s2 at Gata4 weight 2.0, under the ruling row of TRIES.md
    (tries, default TRIES_FILE)."""
    kind, base = parse_form(form)
    ruled = check_gata4(f"{kind}-{base}", gata4_exponent, captain_ruling, tries)
    slot, x, base_sha = BASES[base]
    frozen = check_frozen()
    tree = check_clean()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)  # the frozen builder's own check
    spec = bs.barred_table(dec, rules)[slot]
    if spec["factor"] != FROZEN_FACTOR:
        raise Refused(f"the frozen barred {slot} factor is {spec['factor']}, not {FROZEN_FACTOR}: not built")
    got: dict = {}
    if ruled is None:
        with combined_loader(kind, x, got):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
    else:
        out_p = Path(out)
        if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
            raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
        figs = run_gata4_figures(dec, Path(wt), run_dir, ruled["gata4_exponent"])
        with gata4_loader(ruled["gata4_exponent"], got, figs):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        if got.get("factor_sha256") != figs["factor_sha256"]:
            raise RuntimeError("the factor applied is not the one the line's figures were computed from")
    if got.get("n_loads") != 1 or side["construction"]["factor_exponent"] != 1.0 or \
            side["spec"]["factor"] != FROZEN_FACTOR:
        raise RuntimeError("the build did not apply the combined factor once at w = 1")
    side["inputs"]["factor_note"] = FACTOR_NOTE
    side["combined_build"] = {
        "script": SCRIPT, "script_sha256": bs.sha256(SCRIPT_PATH), "form": f"{kind}-{base}",
        "base": {"name": base, "slot": slot, "x": x, "reference_sha256": base_sha,
                 "rng_streams": side["inputs"]["rng_streams"]},
        "tables": {k: {"file": v["file"], "sha256": v["sha256"], "recorded_in": v["recorded_in"],
                       "pinned_sha256": TABLE_SHA256[k]} for k, v in got["tables"].items()},
        "combination": got["combination"], "where": WHERE, "list_entry": LIST_ENTRY,
        "frozen_files_sha256": frozen, "tree_at_start": tree,
        "identity": f"the {base} base file byte for byte ({base_sha[:12]})" if kind == "identity" else None,
        "depth_step": None,
    }
    if ruled is not None:
        side["combined_build"]["gata4_line"] = {
            "captain_ruling": ruled, "line": G2_LINE, "file": G2_FILE, "gata4_weight": ruled["gata4_exponent"],
            "mab21l2_exponent": 1.0, "factor_rule": "phi_M ** 1 * phi_G ** w, nothing rescaled or clipped after the "
            "product", "r30": G2_R30, "figures": {k: v for k, v in figs.items() if k != "genes"},
            "factor_sha256_applied": got["factor_sha256"], "list_entry": G2_LIST_ENTRY,
        }
    bs.dump_json(side, Path(out).with_suffix(".json"))
    return side


@app.callback()
def main() -> None:
    """The barred file with its Mab21l2 factor combined with the Gata4 table, in one of the fixed forms."""


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
    gata4_exponent: float = typer.Option(None, help="2.0 only, with --form c1-s2 --captain-ruling R22: row R22's "
                                         "second T3 line's file (phi_M * phi_G ** 2), while row R22 of TRIES.md is "
                                         "the captain's order naming it"),
    captain_ruling: str = typer.Option(None, help="R22 only, with --gata4-exponent 2.0 --form c1-s2"),
):
    """Build one fixed form."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        parse_form(form)
        if gata4_exponent is None and captain_ruling is None:
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads)
        else:
            check_gata4(form, gata4_exponent, captain_ruling)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling)
    except FactorMismatch as e:
        bs.say(f"refused (the factor is not the one the line's figures describe): {e}", "red")
        raise typer.Exit(EXIT_FACTOR_MISMATCH) from None
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(2) from None
    chk = s["output"]["check"]
    cb = s["combined_build"]
    bs.say(f"form {cb['form']}: base {cb['base']['name']} (x = {cb['base']['x']}), Gata4 weight "
           f"{cb['combination']['gata4_weight']:.4f}, c = {cb['combination']['c']:.4f}")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


if __name__ == "__main__":
    app()
