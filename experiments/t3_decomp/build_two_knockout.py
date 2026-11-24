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

--composition-from allowed-s1 --captain-ruling R36 (row R36, the captain's named exception of 2026-11-24; off by
default: without both, every form, byte and behaviour above is unchanged and nothing below is read or computed): the
file cD7, T3_c2_s2_cd7.h5ad, which is C2 on the barred s2 (`--form c2-s2`, recorded sha256 47111ac1d07e...) rebuilt by
its own command with one change: its composition block is replaced by submissions.allowed.s1.composition of the
decision file, read by code, so that the five scenario states go to 2^-1.3826521560990648 = 0.3835131214774091 of their
wild-type share (before the composition routine renormalises). Only that block is copied: the same entry's factor
(gata4) and depth step (depth: true, a closed lever) are never taken. C2's factor, row number (max_cells, read by
code), seed, streams and coordinates rule are unchanged, and C2 has no depth step. The replacement is made where the
frozen builder resolves the barred s2's spec (bs.resolve_spec, for this process only), so the frozen validation,
composition routine, draw and assembly run on it as on any spec; at the build the composition routine's result is
checked as G1 checks it, before any cell is drawn.
Refused (exit 2, nothing built): --captain-ruling R36 without the option; any value but allowed-s1; any form but
c2-s2; the option without --captain-ruling R36 or with another ruling; --gata4-exponent beside it; a ledger
(experiments/harness/TRIES.md) without exactly one row R36 whose first cell opens "R36 (the captain's named
exception", says "ruling, no run" and names `--composition-from allowed-s1`, `--form c2-s2` and `--captain-ruling
R36` (the check build_outside_form.captain_order_names makes of an order); a decision file whose sha256 is not the
pinned one (9440ae67...); a copied block whose states differ from rules.yaml's scenario states, or whose log2
magnitude differs from rule D7's m (the decision file's rules.D7.m, 1.3826521560990648), or that is not a central
composition at depletion factor 2 ** -m; a composition routine that does not give the five states at that factor of
their wild-type share. The sidecar's combined_build then also holds `composition_line` (the ruling row, the block and
the one it replaced, the five ratios).

--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22 (row R22's third T3 search line, decision log 2026-11-24 (kk),
on the captain's own yes naming T3_c2_s2_g2.h5ad, ruling row R41; off by default: without the option on c2-s2, every
form, byte and behaviour above is unchanged, combine() included, and nothing below is read or computed): the file
T3_c2_s2_g2.h5ad, C2's command with the gata4_pooled weight g inside C2's size-held product raised from 1 to 2,

  c2-s2 --gata4-exponent 2.0   phi_M ** c * phi_G ** (g c), c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), g = 2.0

population SD over the genes of the two pinned tables as loaded; nothing rescaled or clipped after the product; at
every g the SD of ln factor is phi_M's; g = 0 gives phi_M and g = 1 C2's factor (combine('c2', 1.0)) value for value.
On c2-s2 the option's value is the weight g, not a power: the powers applied are c on transfer_mab21l2_all and g c on
gata4_pooled. Cells, streams, composition, seed, count and board stay C2's (the barred s2's). Refused (exit 2, nothing
built): any other value on c2-s2 (3.0 included); the option without --captain-ruling R22 or with another ruling; a
ledger without exactly one row R22 whose first cell passes build_outside_form.captain_order_names for
"`--form c2-s2 --gata4-exponent 2.0`" (line 2's entry there does not open this line); and, as before, the option on
any form but c1-s2 (line 2) and c2-s2, or either flag alone. The guard, before the build, from released WT E8.75 (den
as build_outside_form.wt_den makes it) and the two pinned tables only: the den-weighted total of the file's factor
lies within 5% of C2's factor's; else exit 6, nothing written. C2's shift from phi_M's total and the file's shift from
phi_M's total are recorded and gate nothing. The factor the build applies must be the very array the guard describes
(its sha256, over WT E8.75's genes); else exit 5, nothing written. The sidecar's combined_build then also holds
`c2_gata4_line` (the ruling row, the line, g, c and the powers applied, the guard's figures and the factor's sha256).

The two target-free gates of row R36's file (they build and write nothing):

  g1 DECISIONS.json --wt WT --composition-from allowed-s1
      before any build: the frozen composition routine (bs.composition), given C2's spec with the block replaced, on
      the reference the build draws from; prints the five states' shares over their wild-type shares before
      renormalising, one per line, and nothing else; exit 0 when each is 0.3835131214774091 within a relative 1e-12,
      1 otherwise, 2 refused, 3 not computed
  g2 DECISIONS.json --wt WT --run-dir RUN --board T3:bcatenin --seed 1000001 --composition-from allowed-s1
      before upload: C2's and cD7's gate bars 0.01 + 0.99 cos^2(dp, 1) on each of row N4's three gene sets, as N4's
      diagnostic computes them (experiments/research_2026_11/diagnostics.py t3-gate: its pseudobulk and cos2_ones; all
      genes, the genes feasible down or up, the top fifth by |dp|), each from the file's own build inputs (the cells
      the build draws, in memory; no built file is read); prints the six bars (C2's three, then cD7's) and nothing
      else; exit 0 when on each set cD7's bar is at most C2's plus 0.01, 1 otherwise, 2 refused, 3 not computed
"""

from __future__ import annotations

import copy
import hashlib
import io
import sys
from contextlib import contextmanager, nullcontext, redirect_stdout
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

# Row R36 (the captain's named exception of 2026-11-24, "Yes, cD7 on C2, my named exception, bound to the read";
# decision log 2026-11-24 (aa) and (bb)): the file cD7, C2 on the barred s2 with its composition block replaced by the
# decision file's submissions.allowed.s1.composition. Off by default; its own code below, apart from row R22's line.
R36_ROW = "R36"
R36_FORM, R36_BASE = "c2-s2", "s2"
R36_OPTION = "`--composition-from allowed-s1`"
R36_NAMES = (R36_OPTION, "`--form c2-s2`", "`--captain-ruling R36`")  # what row R36's first cell must name
R36_OPENER = "the captain's named exception"
R36_FILE = "T3_c2_s2_cd7.h5ad"
C2_S2_SHA256 = "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"  # C2's file (t3two_fb01487)
# --composition-from -> submissions.<branch>.<slot>; only that entry's composition block is copied, never its factor
# (gata4) or its depth step (depth: true; a closed lever): the factor stays C2's and C2 has no depth step
COMPOSITION_FROM = {"allowed-s1": ("allowed", "s1")}
# experiments/t3_decomp/results/gata4/gata4_decisions.json, the run of 2026-10-28 that the frozen rules decided on
R36_DECISION_SHA256 = "9440ae67625608112666d5235e6ceaa5dd9f47a4b0bbd30e16943f32f35c5a8c"
R36_LOG2_MAGNITUDE = 1.3826521560990648  # rule D7's m (the released validation knockout), the block's log2_magnitude
R36_RATIO = 2.0 ** -R36_LOG2_MAGNITUDE  # 0.3835131214774091: each scenario state's share over its wild-type share
R36_RATIO_RTOL = 1e-12  # G1: each of the five ratios within this relative distance of R36_RATIO
R36_N_STATES = 5
R36_G2_MARGIN = 0.01  # G2: on each of N4's three gene sets, cD7's bar at most C2's plus this
N4_SETS = ("all_genes", "feasible_down_or_up", "top20pct_abs_dp")  # row N4's gene sets, in its order
G2_FAILED_LINE = "G2 failed: no upload; C2 stays."  # the launcher's fixed line when G2 does not pass
R36_SOURCE = ("decision file submissions.allowed.s1.composition (--composition-from allowed-s1, row R36), replacing "
              "the barred s2's composition block")
R36_LINE = ("T3 file cD7 (row R36): C2 on the barred s2 (--form c2-s2) rebuilt by its own command with one change, "
            "its composition block replaced by the decision file's submissions.allowed.s1.composition (the five "
            "scenario states at 2^-1.3826521560990648 = 0.3835 of their wild-type share before renormalising); C2's "
            "factor, row number (max_cells, read by code), seed, streams and coordinates rule unchanged")
R36_LIST_ENTRY = "experiments/harness/TRIES.md row R36; decision log 2026-11-24 (aa) and (bb)"
N4_SCRIPT = HERE.parents[1] / "experiments" / "research_2026_11" / "diagnostics.py"  # row N4's t3-gate
EXIT_GATE_FAILED, EXIT_GATE_CRASHED = 1, 3  # g1 and g2: the gate does not pass; it could not be computed

# Row R22's third T3 search line (decision log 2026-11-24 (kk); the captain's own yes naming T3_c2_s2_g2.h5ad, ruling
# row R41): the gata4_pooled weight g inside C2's size-held product on the barred s2. Off by default; its own code
# below, beside combine(), which it leaves unchanged.
C2G_FORM = "c2-s2"
C2G_WEIGHTS = (2.0,)  # --gata4-exponent on c2-s2: the weight g, the single value 2.0
C2G_OPTION = "`--form c2-s2 --gata4-exponent 2.0`"  # what row R22's first cell must name for this line
C2G_FILE = "T3_c2_s2_g2.h5ad"
C2G_MAX_DEN_SHIFT = 0.05  # the guard: the file's den-weighted total within 5% of C2's factor's
C2G_LINE = ("T3 search line 3 (R22): the gata4_pooled weight g in C2's size-held product phi_M ** c * phi_G ** (g c), "
            "c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), on the barred s2; B (g = 0), C2 (g = 1), this file "
            "T3_c2_s2_g2.h5ad (g = 2); the line ends at 2")
C2G_LIST_ENTRY = ("experiments/harness/TRIES.md rows R22 and R41; decision log 2026-11-24 (w), (cc) and (kk)")
C2G_GUARD_RULE = ("the den-weighted total of the file's factor lies within 5% of C2's factor's (den from released WT "
                  "E8.75 as build_outside_form.wt_den makes it, the two pinned tables only); C2's shift and the file's "
                  "shift from phi_M's total are recorded and gate nothing")
EXIT_GUARD_FAILED = 6

Refused = bs.Refused


class FactorMismatch(Refused):
    """The factor the build applies is not the array the line's figures were computed from (exit status 5)."""


class GuardFailed(Refused):
    """Row R22's third line: the file's den-weighted total is not within 5% of C2's (exit status 6, nothing written)."""


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


# ------------------------------------------------------------------ row R36: the file cD7 (off by default)


def ruling_row_names(row: str, opener: str, names, path: Path | None = None) -> bool:
    """Whether the ledger (TRIES_FILE) holds exactly one row `row` whose first cell starts with "<row> (<opener>",
    says "ruling, no run" and names every one of `names` (each given with its backticks): the check
    build_outside_form.captain_order_names makes of a captain's order, for a ruling row that opens otherwise (row R36
    opens "the captain's named exception"). A row that names them only outside its first cell does not pass."""
    p = TRIES_FILE if path is None else Path(path)
    if not p.is_file():
        return False
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    return first.startswith(f"{row} ({opener}") and "ruling, no run" in first and all(n in first for n in names)


def check_composition_from(form, composition_from=None, captain_ruling=None, gata4_exponent=None,
                           tries: Path | None = None) -> dict | None:
    """--composition-from and --captain-ruling R36: None when a command has neither the option nor the ruling R36
    (nothing is read; every such command goes on to check_gata4 as before); else the record, or Refused: the ruling
    R36 without the option, any value but allowed-s1, any form but c2-s2, the option without --captain-ruling R36,
    --gata4-exponent beside it, or a ledger whose row R36 does not pass ruling_row_names for R36_NAMES."""
    if composition_from is None and captain_ruling != R36_ROW:
        return None
    if composition_from is None:
        raise Refused(f"--captain-ruling {R36_ROW} is taken only with --composition-from allowed-s1 on --form "
                      f"{R36_FORM}")
    if not isinstance(composition_from, str) or composition_from not in COMPOSITION_FROM:
        raise Refused(f"--composition-from {composition_from!r} is refused: only allowed-s1")
    if form != R36_FORM:
        raise Refused(f"--composition-from allowed-s1 is taken on --form {R36_FORM} only, not {form!r}")
    if captain_ruling != R36_ROW:
        raise Refused(f"--composition-from allowed-s1 is taken only with --captain-ruling {R36_ROW}, not "
                      f"{captain_ruling!r}")
    if gata4_exponent is not None:
        raise Refused("--composition-from allowed-s1 takes no --gata4-exponent (cD7 is C2 with one change)")
    p = TRIES_FILE if tries is None else Path(tries)
    if not ruling_row_names(R36_ROW, R36_OPENER, R36_NAMES, p):
        raise Refused(f"--composition-from allowed-s1 --captain-ruling {R36_ROW} is refused: "
                      f"experiments/harness/TRIES.md has no single row {R36_ROW} whose first cell is the captain's "
                      f"named exception (\"ruling, no run\") naming {', '.join(R36_NAMES)}")
    return {"row": R36_ROW, "option": R36_OPTION, "composition_from": composition_from,
            "ledger": "experiments/harness/TRIES.md", "ledger_sha256": bs.sha256(p)}


def composition_block(dec: dict, rules: dict, decision: Path, composition_from: str = "allowed-s1") -> tuple[dict,
                                                                                                             dict]:
    """The block --composition-from names, read by code from the decision file (allowed-s1:
    submissions.allowed.s1.composition), and that block only: the entry's factor and depth step are never read.
    Refused unless the file is the pinned one (R36_DECISION_SHA256); the block's depleted and spared states are
    rules.yaml's scenario states, in order (five depleted); its log2 magnitude is rule D7's m, as the decision file
    records it (rules.D7.m) and as pinned here (R36_LOG2_MAGNITUDE); and it is a central composition whose depletion
    factor is 2 ** -magnitude. Returns (a copy of the block, its record)."""
    if not isinstance(composition_from, str) or composition_from not in COMPOSITION_FROM:
        raise Refused(f"--composition-from {composition_from!r} is refused: only allowed-s1")
    branch, slot = COMPOSITION_FROM[composition_from]
    where = f"submissions.{branch}.{slot}.composition"
    got = bs.sha256(Path(decision))
    if got != R36_DECISION_SHA256:
        raise Refused(f"the decision file's sha256 {got[:16]} is not the pinned one ({R36_DECISION_SHA256[:16]}): "
                      "not built")
    try:
        block = dec["submissions"][branch][slot]["composition"]
    except (KeyError, TypeError):
        raise Refused(f"the decision file has no {where}: not built") from None
    if not isinstance(block, dict):
        raise Refused(f"{where} is not a composition block: not built")
    scen = rules["presets"]["gata4"]["scenario"]
    if (block.get("depleted_states") != list(scen["depleted_states"]) or len(scen["depleted_states"]) != R36_N_STATES
            or block.get("spared_states") != list(scen["spared_states"])):
        raise Refused(f"the states of {where} are not rules.yaml's scenario states: not built")
    d7 = (dec.get("rules") or {}).get("D7")
    m = d7.get("m") if isinstance(d7, dict) else None
    l2 = block.get("log2_magnitude")
    if not (type(m) is float and m == R36_LOG2_MAGNITUDE and type(l2) is float and l2 == m):
        raise Refused(f"the log2 magnitude of {where} is not rule D7's m ({R36_LOG2_MAGNITUDE!r}): not built")
    if block.get("kind") != "central" or block.get("depletion_factor") != R36_RATIO:
        raise Refused(f"{where} is not the central composition at depletion factor {R36_RATIO!r}: not built")
    return copy.deepcopy(block), {"source": f"decision file {where}", "decision": str(decision),
                                  "decision_sha256": got, "pinned_sha256": R36_DECISION_SHA256,
                                  "rule_d7_m": m, "taken": "the composition block only (not the entry's factor or "
                                  "depth step)"}


def with_composition(spec: dict, block: dict) -> dict:
    """`spec` with its composition block replaced by a copy of `block`; every other key as it was."""
    return {**spec, "composition": copy.deepcopy(block)}


def composition_ratios(V, f, info: dict, states) -> list[float]:
    """Each of `states`' share in the composition f over its wild-type share V.p, before the routine renormalised:
    f / p on the state divided by f / p on the states the routine left alone (whose one common value is the
    routine's renormalising factor). nan for a state the routine did not find or whose wild-type share is not
    positive, and for every state when the states left alone do not share one value (the routine changed others)."""
    p = np.asarray(V.p, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    nan = [float("nan")] * len(states)
    pos = {t: i for i, t in enumerate(V.labels)}
    found = set(info.get("states_found") or []) - set(info.get("states_missing") or [])
    named = {pos[s] for s in states if s in pos}
    rest = np.array([i for i in range(p.size) if i not in named and p[i] > 0], dtype=int)
    if f.shape != p.shape or not rest.size:
        return nan
    z = f[rest] / p[rest]
    if not (np.all(np.isfinite(z)) and z.min() > 0 and z.max() - z.min() <= 1e-12 * z.max()):
        return nan
    zc = float(z.mean())
    return [float(f[pos[s]] / p[pos[s]] / zc) if s in found and s in pos and p[pos[s]] > 0 else float("nan")
            for s in states]


def ratios_pass(ratios) -> bool:
    """G1's verdict: five ratios, each finite and within R36_RATIO_RTOL of R36_RATIO (0.3835)."""
    r = np.asarray(ratios, dtype=np.float64)
    return (r.shape == (R36_N_STATES,) and bool(np.all(np.isfinite(r)))
            and bool(np.all(np.abs(r - R36_RATIO) <= R36_RATIO_RTOL * R36_RATIO)))


@contextmanager
def r36_composition(block: dict, out: dict):
    """bs.resolve_spec and bs.composition, for this process only, on row R36's line: the frozen resolve_spec, then
    the barred s2's composition block replaced by `block` (any other slot is refused); the frozen composition routine,
    then G1's check of its result (Refused unless the five states are at R36_RATIO of their wild-type share before
    renormalising), before any cell is drawn. `out` receives the block replaced, the ratios and the call counts."""
    orig_resolve, orig_comp = bs.resolve_spec, bs.composition
    states = list(block["depleted_states"])

    def resolve(dec, rules, branch, submission, hashes):
        spec, info = orig_resolve(dec, rules, branch, submission, hashes)
        if (branch, submission) != (BRANCH, R36_BASE):
            raise Refused(f"--composition-from replaces the barred s2's composition only, not {branch} {submission}")
        if set(spec) != {"factor", "composition"} or spec["factor"] != FROZEN_FACTOR:
            raise Refused("the barred s2's spec is not C2's (its factor alone, no depth, stage or coverage step): "
                          "not built")
        out.update(replaced=copy.deepcopy(spec.get("composition")), n_resolves=out.get("n_resolves", 0) + 1)
        return with_composition(spec, block), {**info, "composition_source": R36_SOURCE}

    def comp(V, spec, params):
        f, info = orig_comp(V, spec, params)
        r = composition_ratios(V, f, info, states)
        if not ratios_pass(r):
            raise Refused(f"the composition routine does not give the {R36_N_STATES} scenario states at "
                          f"{R36_RATIO:.4f} of their wild-type share (G1): not built")
        out.update(ratios=dict(zip(states, r, strict=True)), n_compositions=out.get("n_compositions", 0) + 1)
        return f, info

    bs.resolve_spec, bs.composition = resolve, comp
    try:
        yield
    finally:
        bs.resolve_spec, bs.composition = orig_resolve, orig_comp


def run_g1(decision: Path, wt: Path, composition_from: str = "allowed-s1", threads: int = 16) -> list[float]:
    """G1, before any build: the frozen builder's composition routine (bs.composition), given C2's spec (the barred
    s2's, as bs.barred_table gives it) with its composition block replaced as the build replaces it, on the reference
    the build draws from (bs.load_reference on `wt`, the frozen builder's preset); returns the five states' ratios to
    their wild-type shares before renormalising (composition_ratios), in the block's order. Nothing is built or
    written."""
    check_frozen()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    block, _ = composition_block(dec, rules, Path(decision), composition_from)
    spec = with_composition(bs.barred_table(dec, rules)[R36_BASE], block)
    bs.validate_spec(spec, rules)
    ref = bs.load_reference(Path(wt), dec["preset"], rules, None, threads)
    f, info = bs.composition(ref.V, spec, {})
    return composition_ratios(ref.V, f, info, block["depleted_states"])


def _n4():
    """Row N4's diagnostic (experiments/research_2026_11/diagnostics.py), imported by G2 only: its pseudobulk and
    cos2_ones are used as they are, not copied."""
    import importlib.util

    name = "_n4_diagnostics"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, N4_SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        try:
            spec.loader.exec_module(mod)
        except BaseException:
            sys.modules.pop(name, None)
            raise
    return sys.modules[name]


def n4_sets(dp, pb_wt) -> dict[str, np.ndarray]:
    """Row N4's three gene sets for one file's dp, as its t3-gate builds them for each file (N4_SETS): every gene;
    the genes feasible down (pb_wt >= 0.25) or with dp > 0; the fifth of the genes with the largest |dp|."""
    dp = np.asarray(dp, dtype=np.float64)
    feas = np.asarray(pb_wt, dtype=np.float64) >= 0.25
    top = np.argsort(-np.abs(dp))[: max(1, dp.size // 5)]
    sets = {"all_genes": np.arange(dp.size), "feasible_down_or_up": np.flatnonzero(feas | (dp > 0)),
            "top20pct_abs_dp": top}
    return {s: sets[s] for s in N4_SETS}


def n4_bars(dp, pb_wt) -> dict[str, float]:
    """Row N4's gate bar 0.01 + 0.99 cos^2(dp, 1) on each of its three gene sets (n4_sets), cos2_ones its own."""
    dp = np.asarray(dp, dtype=np.float64)
    return {s: 0.01 + 0.99 * _n4().cos2_ones(dp[idx]) for s, idx in n4_sets(dp, pb_wt).items()}


def g2_passes(bars_c2: dict, bars_cd7: dict) -> bool:
    """G2's verdict: on each of N4's three gene sets both bars are finite and cD7's is at most C2's plus
    R36_G2_MARGIN."""
    return all(s in bars_c2 and s in bars_cd7 and bool(np.isfinite(bars_c2[s]) and np.isfinite(bars_cd7[s]))
               and bars_cd7[s] <= bars_c2[s] + R36_G2_MARGIN for s in N4_SETS)


def run_g2(decision: Path, wt: Path, run_dir: Path | None, board, n_cells: int, seed: int = 1_000_001,
           composition_from: str = "allowed-s1", threads: int = 16) -> dict:
    """G2, before upload: C2's and cD7's gate bars on N4's three gene sets (n4_bars), as row N4's diagnostic computes
    them: dp is the pseudobulk of the drawn cells over the board's panel minus that of `wt` as released
    (vec.prep.load_stage, as the scorer loads the reference), both by N4's pseudobulk. Each file's cells come from its
    own build inputs, through the build's own code up to the draw: bs.preflight on the barred s2 (cD7's under
    r36_composition), bs.load_reference, the run records check, bs.inputs_for with C2's factor (combined_loader, c2 at
    x = 1) and bs.construct at the seed and the stream "submission:barred:s2". Nothing is assembled, written or read
    from a built file. Returns both files' three bars and the verdict (g2_passes)."""
    from vec.prep import load_stage

    n4 = _n4()
    check_frozen()
    rules = bs.load_rules()
    dec = bs.load_decision(Path(decision), rules)
    bs.check_hashes(dec, allow_mismatch=False)
    block, _ = composition_block(dec, rules, Path(decision), composition_from)
    kind, base = parse_form(R36_FORM)
    slot, x, _ = BASES[base]
    ref, dps, hook = None, {}, {}
    for name in ("c2", "cd7"):
        got: dict = {}
        hooked = r36_composition(block, hook) if name == "cd7" else nullcontext()
        with hooked, combined_loader(kind, x, got):
            pf = bs.preflight(Path(decision), BRANCH, slot, board, n_cells)
            if ref is None:
                ref = bs.load_reference(Path(wt), pf["preset"], rules, None, threads)
                bs.check_inputs_against_run(bs.run_records(run_dir, pf["preset"]), pf["preset"],
                                            ref.record["sha256"], {})
                panel = list(board.genes)
                gpos = {g: i for i, g in enumerate(ref.genes)}
                missing = [g for g in panel if g not in gpos]
                if missing:
                    raise Refused(f"{len(missing)} board genes are not in the WT reference")
                cols = [gpos[g] for g in panel]
                wt_rel = load_stage(Path(wt))
                wpos = {g: i for i, g in enumerate(wt_rel.var_names)}
                pb_wt = n4.pseudobulk(wt_rel.X.tocsr()[:, [wpos[g] for g in panel]])
            kw, _ = bs.inputs_for(pf["spec"], pf["dec"], ref, run_dir)
            built = bs.construct(ref, n_cells, bs.t3d.stable_rng(seed, f"submission:{BRANCH}:{slot}"), **kw)
        if got.get("n_loads") != 1 or kw["w"] != 1.0:
            raise RuntimeError("the draw did not apply C2's combined factor once at w = 1")
        dps[name] = n4.pseudobulk(built.X[:, cols]) - pb_wt
        del built
    if hook.get("n_resolves") != 1 or hook.get("n_compositions") != 1:
        raise RuntimeError("cD7's composition was not replaced and checked exactly once")
    bars = {name: n4_bars(dp, pb_wt) for name, dp in dps.items()}
    return {"c2": bars["c2"], "cd7": bars["cd7"], "margin": R36_G2_MARGIN, "sets": list(N4_SETS),
            "passes": g2_passes(bars["c2"], bars["cd7"])}


# ------------------------------------------------------------- row R22's third T3 line: T3_c2_s2_g2 (off by default)


def check_c2_gata4(form, gata4_exponent=None, captain_ruling=None, tries: Path | None = None) -> dict | None:
    """--gata4-exponent on --form c2-s2 (row R22's third T3 line): None unless the form is c2-s2 and the option is
    given (nothing is read; every other command, the ruling alone on c2-s2 included, goes on to check_gata4 as
    before); else the line's record, or Refused: a ruling other than R22, a value other than 2.0, the option without
    the ruling, or a ledger whose row R22 does not pass build_outside_form.captain_order_names for C2G_OPTION (line
    2's entry there, `--gata4-exponent 2.0` on `--form c1-s2`, does not open this line)."""
    if form != C2G_FORM or gata4_exponent is None:
        return None
    if captain_ruling is not None and captain_ruling != G2_ROW:
        raise Refused(f"--captain-ruling {captain_ruling!r} is refused: only row {G2_ROW} orders a file of this "
                      "builder")
    try:
        g = float(gata4_exponent)
    except (TypeError, ValueError):
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0") from None
    if isinstance(gata4_exponent, bool) or g not in C2G_WEIGHTS:
        raise Refused(f"--gata4-exponent {gata4_exponent!r} is refused: only 2.0")
    if captain_ruling is None:
        raise Refused(f"--gata4-exponent 2.0 on --form {C2G_FORM} is taken only with --captain-ruling {G2_ROW}")
    p = TRIES_FILE if tries is None else Path(tries)
    if not _outside_form().captain_order_names(G2_ROW, C2G_OPTION, p):
        raise Refused(f"--form {C2G_FORM} --gata4-exponent 2.0 --captain-ruling {G2_ROW} is refused: "
                      f"experiments/harness/TRIES.md has no single row {G2_ROW} whose first cell is the captain's "
                      f"order (\"ruling, no run\") naming {C2G_OPTION} and `--captain-ruling {G2_ROW}`")
    return {"row": G2_ROW, "option": C2G_OPTION, "gata4_weight": g, "ledger": "experiments/harness/TRIES.md",
            "ledger_sha256": bs.sha256(p)}


def c2_gata4_weighted(phi_m, phi_g, g: float) -> tuple[np.ndarray, dict]:
    """C2's size-held product with the gata4_pooled weight g, from the two tables as loaded, gene by gene:
    phi_M ** c * phi_G ** (g c), c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD over the genes; nothing
    is rescaled or clipped after the product. g = 0 gives phi_M and g = 1 combine('c2', 1.0, ...)'s array, value for
    value (the same operations in the same order); at every g the SD of ln factor is phi_M's. Returns (phi, record);
    the record names g (the option's value on c2-s2: a weight, not a power) and the powers applied, c and g c."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    phi_g = np.asarray(phi_g, dtype=np.float64)
    if phi_m.shape != phi_g.shape or phi_m.ndim != 1:
        raise Refused("the two factor tables do not cover the same genes")
    for nm, p in ((TABLE_M, phi_m), (TABLE_G, phi_g)):
        if not np.all(np.isfinite(p) & (p > 0)):
            raise Refused(f"{nm}: factor values must be finite and positive")
    g = float(g)
    if not np.isfinite(g):
        raise Refused(f"the gata4_pooled weight {g!r} is not finite")
    lm, lg = np.log(phi_m), np.log(phi_g)
    sd_base, sd_sum = _sd(lm), _sd(lm + g * lg)
    if not sd_sum > 0:
        raise Refused("SD(ln phi_M + g ln phi_G) is 0: c is undefined")
    c = sd_base / sd_sum
    phi = (phi_m ** c) * (phi_g ** (g * c))
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused("the combined factor is not finite and positive")
    sd_m, sd_g = _sd(lm), _sd(lg)
    rec = {"kind": "c2", "x": 1.0, "g": g, "gata4_weight": g * c, "mab21l2_exponent": float(c), "c": float(c),
           "c_rule": "SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD over the genes; on c2-s2 the option's "
                     "value is the weight g, not a power: the powers applied are c (transfer_mab21l2_all) and g c "
                     "(gata4_pooled)",
           "sd_ln_phi_m": sd_m, "sd_ln_phi_g": sd_g, "sd_ln_base": sd_base, "sd_ln_base_plus_wg": sd_sum,
           "corr_ln_m_ln_g": float(np.corrcoef(lm, lg)[0, 1]) if sd_m > 0 and sd_g > 0 else None,
           "base_factor": _shift(phi_m), "factor": _shift(phi),
           "den_weighted_mean": "recorded in combined_build.c2_gata4_line.guard",
           "min": float(phi.min()), "max": float(phi.max()), "n_up": int((phi > 1).sum()),
           "n_down": int((phi < 1).sum())}
    return phi, rec


def c2g_guard(phi_m, phi_g, g: float, den) -> dict:
    """The line's guard, a pure function of the two tables and den: the den-weighted totals (den . phi / sum den) of
    phi_M (B's factor), of C2's factor (combine('c2', 1.0), as C2's build applies it) and of the file's
    (c2_gata4_weighted at g); passes when the file's total lies within C2G_MAX_DEN_SHIFT of C2's. C2's shift and the
    file's shift from phi_M's total are recorded and gate nothing. Also the sha256 of the factor it describes."""
    phi, rec = c2_gata4_weighted(phi_m, phi_g, g)
    c2, _ = combine("c2", 1.0, phi_m, phi_g)
    phi_m = np.asarray(phi_m, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    if den.shape != phi_m.shape or not np.all(np.isfinite(den) & (den >= 0)) or not den.sum() > 0:
        raise Refused("den does not cover the tables' genes with finite, non-negative weights and a positive sum")
    s = float(den.sum())
    t_m, t_c2, t_f = (float(den @ f / s) for f in (phi_m, c2, phi))
    shift = t_f / t_c2 - 1.0
    return {"gata4_weight_g": float(g), "c": rec["c"],
            "powers_applied": {TABLE_M: rec["mab21l2_exponent"], TABLE_G: rec["gata4_weight"]},
            "den_total_phi_m": t_m, "den_total_c2": t_c2, "den_total_factor": t_f, "den_shift_vs_c2": shift,
            "max_den_shift": C2G_MAX_DEN_SHIFT, "passed": bool(abs(shift) <= C2G_MAX_DEN_SHIFT),
            "den_shift_c2_vs_phi_m": t_c2 / t_m - 1.0, "den_shift_factor_vs_phi_m": t_f / t_m - 1.0,
            "vs_phi_m": "recorded; gates nothing", "rule": C2G_GUARD_RULE, "factor_sha256": factor_sha256(phi)}


def run_c2g_guard(dec: dict, wt: Path, run_dir: Path | None, g: float) -> dict:
    """The guard from released WT E8.75 and the two pinned tables only, before any build: den and genes by
    build_outside_form.wt_den, both tables through the frozen loader in that gene order (checked against
    TABLE_SHA256). Returns the guard's record whatever it says; the caller refuses on it."""
    genes, den = _outside_form().wt_den(Path(wt))
    phi_m, rec_m = bs.load_factor(run_dir, TABLE_M, genes, dec)
    phi_g, rec_g = bs.load_factor(run_dir, TABLE_G, genes, dec)
    for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
        if rec["sha256"] != TABLE_SHA256[nm]:
            raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                          f"({TABLE_SHA256[nm][:16]}): not built")
    den = np.asarray(den, dtype=np.float64)
    out = c2g_guard(phi_m, phi_g, g, den)
    out |= {"wt": str(wt), "wt_sha256": bs.sha256(Path(wt)), "den_rule": G2_DEN_RULE, "den": den.tolist(),
            "den_sha256": factor_sha256(den), "n_genes": len(genes),
            "genes_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
            "tables_sha256": {TABLE_M: rec_m["sha256"], TABLE_G: rec_g["sha256"]}, "genes": list(genes)}
    return out


@contextmanager
def c2g_loader(g: float, out: dict, guard: dict):
    """bs.load_factor, for this process only, on row R22's third line: the frozen loader for phi_M and phi_G (checked
    against TABLE_SHA256 too), then c2_gata4_weighted at g, the hook C2's factor goes through. FactorMismatch (exit 5)
    unless the build's genes are the guard's and the factor's sha256 is the one the guard describes."""
    original = bs.load_factor

    def load(run_dir, name, genes, dec):
        if name != TABLE_M:
            raise Refused(f"the barred spec asked for the factor {name!r}, not {TABLE_M!r}: not built")
        if list(genes) != guard["genes"]:
            raise FactorMismatch("the build's genes are not those the guard was computed over: not built")
        phi_m, rec_m = original(run_dir, TABLE_M, genes, dec)
        phi_g, rec_g = original(run_dir, TABLE_G, genes, dec)
        for nm, rec in ((TABLE_M, rec_m), (TABLE_G, rec_g)):
            if rec["sha256"] != TABLE_SHA256[nm]:
                raise Refused(f"{rec['file']}: sha256 {rec['sha256'][:16]} differs from the pinned table "
                              f"({TABLE_SHA256[nm][:16]}): not built")
        phi, crec = c2_gata4_weighted(phi_m, phi_g, g)
        sha = factor_sha256(phi)
        if sha != guard["factor_sha256"]:
            raise FactorMismatch(f"the factor the build would apply ({sha[:16]}) is not the one the guard describes "
                                 f"({guard['factor_sha256'][:16]}): not built")
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
               captain_ruling=None, tries: Path | None = None, composition_from=None) -> dict:
    """The named form through bs.build_submission; returns the sidecar (written next to the file, the frozen
    builder's plus `combined_build`). gata4_exponent 2.0 with captain_ruling 'R22' (default None, today's paths
    unchanged) builds row R22's second T3 line's file, c1-s2 at Gata4 weight 2.0, under the ruling row of TRIES.md
    (tries, default TRIES_FILE). composition_from 'allowed-s1' with captain_ruling 'R36' (default None) builds row
    R36's file cD7, c2-s2 with the decision file's allowed s1 composition block, under the ruling row of TRIES.md.
    gata4_exponent 2.0 with captain_ruling 'R22' on c2-s2 builds row R22's third T3 line's file, C2 at gata4_pooled
    weight g = 2.0, after its guard, under the ruling row of TRIES.md."""
    kind, base = parse_form(form)
    cd7 = check_composition_from(f"{kind}-{base}", composition_from, captain_ruling, gata4_exponent, tries)
    c2g = None if cd7 is not None else check_c2_gata4(f"{kind}-{base}", gata4_exponent, captain_ruling, tries)
    ruled = None if cd7 is not None or c2g is not None else check_gata4(f"{kind}-{base}", gata4_exponent,
                                                                         captain_ruling, tries)
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
    if cd7 is not None:
        block, brec = composition_block(dec, rules, Path(decision), composition_from)
        hook: dict = {}
        with combined_loader(kind, x, got), r36_composition(block, hook):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        if hook.get("n_resolves") != 1 or hook.get("n_compositions") != 1 or side["spec"]["composition"] != block:
            raise RuntimeError("the build did not replace and check the composition block exactly once")
        if set(side["spec"]) != {"factor", "composition"} or side["construction"]["depth_ratios"] is not None:
            raise RuntimeError("the build took more than the composition block")
    elif c2g is not None:
        out_p = Path(out)
        if not overwrite and (out_p.exists() or out_p.with_suffix(".json").exists()):
            raise Refused(f"{out_p} or its sidecar exists; pass --overwrite")
        guard = run_c2g_guard(dec, Path(wt), run_dir, c2g["gata4_weight"])
        if not guard["passed"]:
            raise GuardFailed(f"the file's den-weighted total is {guard['den_shift_vs_c2']:+.2%} against C2's "
                              f"factor's, beyond {C2G_MAX_DEN_SHIFT:.0%}: not built")
        with c2g_loader(c2g["gata4_weight"], got, guard):
            side = bs.build_submission(Path(decision), BRANCH, slot, Path(wt), board, n_cells, Path(out),
                                       run_dir=run_dir, seed=seed, overwrite=overwrite, threads=threads)
        if got.get("factor_sha256") != guard["factor_sha256"]:
            raise RuntimeError("the factor applied is not the one the guard describes")
    elif ruled is None:
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
    if cd7 is not None:
        side["combined_build"]["composition_line"] = {
            "captain_ruling": cd7, "line": R36_LINE, "file": R36_FILE, "composition_from": composition_from,
            "decision_file": brec, "composition": block, "replaced": hook["replaced"],
            "log2_magnitude": R36_LOG2_MAGNITUDE, "ratio": R36_RATIO,
            "ratios_before_renormalising": hook["ratios"], "c2_reference_sha256": C2_S2_SHA256,
            "unchanged": "C2's factor (combined_build.combination), row number (max_cells, read by code), seed, "
                         "streams and coordinates rule",
            "list_entry": R36_LIST_ENTRY,
        }
    if c2g is not None:
        crec = got["combination"]
        side["combined_build"]["c2_gata4_line"] = {
            "captain_ruling": c2g, "line": C2G_LINE, "file": C2G_FILE, "gata4_weight_g": c2g["gata4_weight"],
            "c": crec["c"], "powers_applied": {TABLE_M: crec["mab21l2_exponent"], TABLE_G: crec["gata4_weight"]},
            "option_value": "on --form c2-s2 the value of --gata4-exponent is the gata4_pooled weight g, not a power",
            "factor_rule": "phi_M ** c * phi_G ** (g c), c = SD(ln phi_M) / SD(ln phi_M + g ln phi_G), population SD "
                           "over the genes; nothing rescaled or clipped after the product",
            "guard": {k: v for k, v in guard.items() if k != "genes"},
            "factor_sha256_applied": got["factor_sha256"], "c2_reference_sha256": C2_S2_SHA256,
            "list_entry": C2G_LIST_ENTRY,
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
    gata4_exponent: float = typer.Option(None, help="2.0 only, with --captain-ruling R22: on --form c1-s2 row R22's "
                                         "second T3 line's file (phi_M * phi_G ** 2); on --form c2-s2 its third "
                                         "line's file, the gata4_pooled weight g = 2 in C2's size-held product (a "
                                         "weight, not a power); each while row R22 of TRIES.md is the captain's order "
                                         "naming it"),
    captain_ruling: str = typer.Option(None, help="R22, with --gata4-exponent 2.0 --form c1-s2 or c2-s2; or R36, "
                                       "with --composition-from allowed-s1 --form c2-s2; nothing else"),
    composition_from: str = typer.Option(None, help="allowed-s1 only, with --form c2-s2 --captain-ruling R36: row "
                                         "R36's file cD7 (C2 with the decision file's submissions.allowed.s1."
                                         "composition), while row R36 of TRIES.md is the captain's ruling naming it"),
):
    """Build one fixed form."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        parse_form(form)
        if composition_from is not None or captain_ruling == R36_ROW:
            check_composition_from(form, composition_from, captain_ruling, gata4_exponent)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, captain_ruling=captain_ruling,
                           composition_from=composition_from)
        elif form == C2G_FORM and gata4_exponent is not None:
            check_c2_gata4(form, gata4_exponent, captain_ruling)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling)
        elif gata4_exponent is None and captain_ruling is None:
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads)
        else:
            check_gata4(form, gata4_exponent, captain_ruling)
            b = bs.board_spec(board, boards_dir)
            s = build_form(decision, form, wt, b, bs.resolve_n_cells(n_cells, b), out, run_dir=run_dir, seed=seed,
                           overwrite=overwrite, threads=threads, gata4_exponent=gata4_exponent,
                           captain_ruling=captain_ruling)
    except GuardFailed as e:
        bs.say(f"refused (the guard: the file's den-weighted total against C2's): {e}", "red")
        raise typer.Exit(EXIT_GUARD_FAILED) from None
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
    if cb.get("c2_gata4_line"):
        bs.say(f"row {G2_ROW}'s third T3 line: gata4_pooled weight g = {cb['c2_gata4_line']['gata4_weight_g']:g} in "
               f"C2's size-held product (a weight, not a power); the guard passed")
    if cb.get("composition_line"):
        bs.say(f"composition from {cb['composition_line']['composition_from']} (row {R36_ROW}): the "
               f"{R36_N_STATES} scenario states at {R36_RATIO:.4f} of their wild-type share before renormalising")
    for w in s["warnings"]:
        bs.say(f"  warning: {w}", "yellow")
    for e in chk["errors"]:
        bs.say(f"  check error: {e}", "red")
    bs.say(f"{'OK' if chk['ok'] else 'FAILED the format check'}: {s['output']['h5ad']}",
           "green" if chk["ok"] else "red")
    raise typer.Exit(0 if chk["ok"] else 1)


@app.command()
def g1(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the build"),
    composition_from: str = typer.Option(..., help="allowed-s1 only"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """G1, row R36's gate before any build: prints the five scenario states' shares over their wild-type shares
    before renormalising, one per line in the block's order (Python's shortest round-trip form), and nothing else;
    exit 0 when each is 0.3835131214774091 within a relative 1e-12, 1 otherwise, 2 refused, 3 not computed."""
    try:
        with redirect_stdout(io.StringIO()):  # anything the loading code prints is dropped: the ratios only
            ratios = run_g1(decision, wt, composition_from, threads)
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    except Exception as e:  # the type only: a message or traceback could carry a figure
        bs.say(f"failed: {type(e).__name__}", "red")
        raise typer.Exit(EXIT_GATE_CRASHED) from None
    for r in ratios:
        print(repr(float(r)))
    raise typer.Exit(0 if ratios_pass(ratios) else EXIT_GATE_FAILED)


@app.command()
def g2(
    decision: Path = typer.Argument(..., help="decision file (experiments/t3_decomp/results/gata4/"
                                    "gata4_decisions.json)"),
    wt: Path = typer.Option(..., help="WT reference (WT E8.75), as for the build"),
    board: str = typer.Option(..., help="T3 board key in boards/"),
    run_dir: Path = typer.Option(None, help="run directory holding factors/"),
    seed: int = typer.Option(1_000_001, help="1000001, the list's seed; any other is refused"),
    n_cells: str = typer.Option("max", help="'max', the build's count: the board's max_cells, read in code, never "
                                "printed; any other is refused"),
    boards_dir: Path = typer.Option(None, help="boards directory (default: the repository's boards/)"),
    composition_from: str = typer.Option(..., help="allowed-s1 only"),
    threads: int = typer.Option(16, help="BLAS threads while preparing the reference (the frozen default)"),
):
    """G2, row R36's gate before upload: prints the six gate bars, one per line, and nothing else: C2's on N4's three
    gene sets (all genes, feasible down or up, the top fifth by |dp|), then cD7's on the same three; exit 0 when on
    each set cD7's bar is at most C2's plus 0.01, 1 otherwise, 2 refused, 3 not computed."""
    try:
        if seed != LIST_SEED or n_cells != LIST_N_CELLS:
            raise Refused(LIST_INPUTS_REFUSAL)
        if not isinstance(composition_from, str) or composition_from not in COMPOSITION_FROM:
            raise Refused(f"--composition-from {composition_from!r} is refused: only allowed-s1")
        b = bs.board_spec(board, boards_dir)
        with redirect_stdout(io.StringIO()):  # anything the loading code prints is dropped: the six bars only
            res = run_g2(decision, wt, run_dir, b, bs.resolve_n_cells(n_cells, b), seed, composition_from, threads)
    except Refused as e:
        bs.say(f"refused: {e}", "red")
        raise typer.Exit(EXIT_REFUSED) from None
    except Exception as e:  # the type only: a message or traceback could carry a count of the test build
        bs.say(f"failed: {type(e).__name__}", "red")
        raise typer.Exit(EXIT_GATE_CRASHED) from None
    for name in ("c2", "cd7"):
        for s in N4_SETS:
            print(repr(float(res[name][s])))
    raise typer.Exit(0 if res["passes"] else EXIT_GATE_FAILED)


if __name__ == "__main__":
    app()
