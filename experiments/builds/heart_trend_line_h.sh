#!/usr/bin/env bash
# The heart's trend-form line (row R55 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling, no
# run: his named exception on the heart's scope and size for OTK125, OTX125, OTW125 and OTP125, bound to their screen,
# row H35; the HEART PORTFOLIO design of the step-change round as its cold challenge fixes it, OKS125 struck): OTS125's
# command (file D at strength 1.25, row R49's) with one change, --trend-form <F> under --captain-ruling R55:
#   category  OTK125  T2heart_organ_ts125_tk.h5ad      contrast  OTX125  T2heart_organ_ts125_tx.h5ad
#   whole     OTW125  T2heart_organ_ts125_tw.h5ad      portion   OTP125  T2heart_organ_ts125_tp.h5ad
# each then halved by experiments/builds/halve.py (seed 1000001) as OTS125 is. A form is built only when row H35's
# outcome cell, read by code, says KEPT for it: the cell must quote the read's outcome line, "OUTCOME H35: category
# <W>; contrast <W>; whole <W>; portion <W>." (W one of KEPT, DROP and IDLE); a VOID or STOP outcome, a missing or
# malformed row H35 or an outcome cell without that line builds nothing. Test-board builds on a DGX campaign checkout
# (never the reference checkout). Run from the root of the checkout or detached worktree it is called in, detached, and
# wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_trend_line_h.sh <COMMIT> <RUN>
# The row is fixed (R55; the builder pins it), so there is no ROW argument. RUN names this run's own scratch folder
# ~/v/scratch/hearth/<RUN> and its output folder runs/hearth_<commit7>_<RUN>. Neither may exist: a run folder is never
# reused and nothing is ever deleted here, so a second run takes a new name. In order, each step writing
# ~/v/scratch/hearth/<RUN>/<STEP>.exit (the word running first, once its checks pass, then its exit code):
#   gate    HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks; TRIES.md holds
#           exactly one row H35 of eleven cells whose first cell is the screen's fixed text and whose outcome cell
#           allows at least one form (the forms it allows go to the scratch folder's forms file); and exactly one row
#           R55 whose first cell is the captain's order ("ruling, no run") naming `--entry organ_scope`,
#           `--trend-strength 1.25`, `--trend-form <F>` for each allowed form and `--captain-ruling R55`, each as a
#           whole backticked string (the builder's own check); then the screen's run commit, read by code from its
#           records as committed at COMMIT (experiments/harness/results/h35_heart_trend_forms: G0's g0.json and the
#           read's reading.json must name one commit, G0's from a clean tree, and the read must keep exactly the forms
#           row H35's outcome cell keeps), must be in COMMIT's history, and no file the screen read (SCREEN_FILES: src,
#           configs, the builder, halve.py and the screen's scripts) may differ between it and COMMIT. So the files
#           are built from the code the screen read: build on the screen's run commit with only the ledger's and the
#           records' commits on top, or have the change reviewed and screened again. The run commit goes to the
#           scratch folder's screened file
#   tests   the trend-form line's, its screen's and the earlier heart lines' targeted tests (synthetic stages only),
#           8 BLAS threads
#   D OTS05 OTS10 OTS075 OTS125 OTS15   the identity steps, OTS15's included: each file by its committed command,
#           then its half, and both must give their recorded sha256 (experiments/builds/results/
#           heartg_b2bbcf5_g1.SHA256SUMS, the last build of all six), or everything stops
#   <tag>a <tag>b   each allowed form's file, twice, into separate folders, each then halved by halve.py
#   done    one sha256 list over every built file; each form's two builds compared byte for byte, file and half
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or
# sidecar record. A build that exits non-zero means no file, and the reason is not read here. Nothing here uploads
# anything; a file goes in only on the captain's word at the time.
# Exit codes: 90 the arguments, 91 the run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 the
# rules pin, 92 row H35 missing or malformed or its outcome cell without the read's outcome line, 87 the outcome is
# VOID, 86 the outcome is STOP, 95 row R55 does not name an allowed form's file, 84 the screen's records do not name
# one clean run commit in COMMIT's history that keeps the same forms, 85 a file the screen read differs between that
# commit and COMMIT, 93 the output folder exists, 10 the
# tests, 11 a build failed, 12 halve.py, 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 20 to 29 and
# 32 33 an identity step is not its recorded bytes (file, half: D 20 21, OTS05 22 23, OTS10 24 25, OTS075 26 27, OTS125
# 28 29, OTS15 32 33), 30 the hash list, 31 a form's two builds differ.
set -u
set -o pipefail
COMMIT=${1:?usage: heart_trend_line_h.sh COMMIT RUN}
RUN=${2:?usage: heart_trend_line_h.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
[ "$#" = 2 ] || { echo "refused: COMMIT and RUN only (the row, R55, is fixed)" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/hearth/$RUN
mkdir -p "$HOME/v/scratch/hearth"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
D_SHA=76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b            # heartg g1: T2heart_organ.h5ad
D_HALF_SHA=1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00       # heartg g1: its half
OTS05_SHA=c3b5fc9c7b47d0d0aafc4e5876bac7f06988cb2d628a0030b4b643806da8bceb        # heartg g1: OTS05's file
OTS05_HALF_SHA=f918e76731aac18eda3b46c8a35c18f26ae00a1a6a1db36d33e7515704cc43ea   # heartg g1: OTS05's half
OTS10_SHA=efaaabe0b744d82dc40135ed9269a6baf8938c9cdd83f85fb3917c0b102e7a2b        # heartg g1: OTS10's file
OTS10_HALF_SHA=ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9   # heartg g1: OTS10's half
OTS075_SHA=109f7f26966c4a3cd0cd18807032ea8eb8f7786df833cd7b3064a21298b1cdbb       # heartg g1: OTS075's file
OTS075_HALF_SHA=1a66b5a5a9d5388d806f548d5c51570345f02e340a6ef2d4ce3490d1e2673f5d  # heartg g1: OTS075's half
OTS125_SHA=79d0b813db758ce88bda22ed1f8c8b58a560214360bfbfc2795c7eee372419da       # heartg g1: OTS125's file
OTS125_HALF_SHA=adc3243665949f9802056216d86443e2ea0d45d07534013f99e26afb7e57ade0  # heartg g1: OTS125's half
OTS15_SHA=ade351c1627d0eb5bdbc43cb613cf1b914053e19e4efca0dd943f939b1035d08        # heartg g1: OTS15's file
OTS15_HALF_SHA=f818931834d08a6906e010e8af9f077ba15be779aba8965cbaa880a9bdf1196c   # heartg g1: OTS15's half
D_ARGS="--entry organ_scope --expression copy --z-mode raw"
OTS05_ARGS="--entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29"
OTS10_ARGS="--entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29"
OTS075_ARGS="--entry organ_scope --trend-strength 0.75 --z-mode raw --captain-ruling R49"
OTS125_ARGS="--entry organ_scope --trend-strength 1.25 --z-mode raw --captain-ruling R49"
OTS15_ARGS="--entry organ_scope --trend-strength 1.5 --z-mode raw --captain-ruling R52"
TK_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55"
TX_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form contrast --z-mode raw --captain-ruling R55"
TW_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form whole --z-mode raw --captain-ruling R55"
TP_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form portion --z-mode raw --captain-ruling R55"
# The screen's committed records, and the files the screen read (its scripts as they import each other, the builder
# whose context G0 takes, halve.py, the size procedure, every vec module and every config)
SCREEN_RECORDS=experiments/harness/results/h35_heart_trend_forms
SCREEN_FILES="src configs experiments/builds/build_test.py experiments/builds/halve.py"
SCREEN_FILES="$SCREEN_FILES experiments/coords/size_procedure.py experiments/s2_hedge/payoff.py"
SCREEN_FILES="$SCREEN_FILES experiments/harness/h25_heart_forms.py experiments/harness/h26_heart_scope_forms.py"
SCREEN_FILES="$SCREEN_FILES experiments/harness/h27_heart_forecast_size.py"
SCREEN_FILES="$SCREEN_FILES experiments/harness/h28_heart_trend_strength.py experiments/harness/h29b_heart_organ.py"
SCREEN_FILES="$SCREEN_FILES experiments/harness/h34_heart_organ_trend.py experiments/harness/h35_heart_trend_forms.py"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
FORMS=$(python3 -I -B -c '
import re
import sys
FIRST = "H35 (the counted screen of row R55\x27s trend forms: ten tries; G0 first, target-free, not a try)"
NAMES = ("category", "contrast", "whole", "portion")
MARK = "OUTCOME H35:"
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| H35 ")]
line = rows[0].rstrip() if len(rows) == 1 else ""
cells = [c.strip() for c in line[1:-1].split("|")] if line.endswith("|") else []
if len(cells) != 11 or cells[0] != FIRST:
    sys.exit(3)
cell = cells[10]
if cell.count(MARK) != 1:
    sys.exit(4)
if MARK + " VOID" in cell:
    sys.exit(5)
words = "(KEPT|DROP|IDLE)"
m = re.search(MARK + " (STOP: )?" + "; ".join(n + " " + words for n in NAMES) + r"\.", cell)
if m is None:
    sys.exit(4)
kept = [n for n, w in zip(NAMES, m.groups()[1:]) if w == "KEPT"]
if bool(m.group(1)) == bool(kept):
    sys.exit(4)
if not kept:
    sys.exit(6)
print(" ".join(kept))
' experiments/harness/TRIES.md)
code=$?
case "$code" in
  0) ;;
  5) echo "refused: row H35's outcome is VOID; nothing is built" >&2; finish 87 ;;
  6) echo "refused: row H35's outcome is STOP; nothing is built" >&2; finish 86 ;;
  *) echo "refused: row H35 is missing or malformed, or its outcome cell holds no outcome line of the read" >&2
     finish 92 ;;
esac
for form in $FORMS; do
  python3 -I -B -c 'import re
import sys
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R55 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
named = set(re.findall(r"`([^`]*)`", first))
need = ("--entry organ_scope", "--trend-strength 1.25", "--trend-form " + sys.argv[2], "--captain-ruling R55")
sys.exit(0 if first.startswith("R55 (the captain\x27s order") and "ruling, no run" in first
         and all(s in named for s in need) else 1)' experiments/harness/TRIES.md "$form" \
    || { echo "refused: TRIES.md does not hold exactly one row R55 that is the captain's order naming $form" >&2
         finish 95; }
done
SCREENED=$(python3 -I -B -c '
import json
import re
import sys
def load(path):
    try:
        with open(path, encoding="utf-8") as f:
            rec = json.load(f)
    except (OSError, ValueError):
        sys.exit(1)
    if not isinstance(rec, dict):
        sys.exit(1)
    return rec
g0, read = load(sys.argv[1]), load(sys.argv[2])
commit = g0.get("commit")
prov = g0.get("provenance") if isinstance(g0.get("provenance"), dict) else {}
words = read.get("words") if isinstance(read.get("words"), dict) else {}
kept = " ".join(n for n in ("category", "contrast", "whole", "portion") if words.get(n) == "KEPT")
ok = (isinstance(commit, str) and re.fullmatch("[0-9a-f]{40}", commit) is not None and prov.get("commit") == commit
      and prov.get("tree_clean") is True and read.get("commit") == commit and read.get("void") is False
      and read.get("outcome") == "KEPT" and kept == sys.argv[3])
if not ok:
    sys.exit(1)
print(commit)
' "$SCREEN_RECORDS/g0.json" "$SCREEN_RECORDS/reading.json" "$FORMS") \
  || { echo "refused: the screen's records do not name one clean run commit that keeps these forms" >&2; finish 84; }
git merge-base --is-ancestor "$SCREENED" "$COMMIT" > /dev/null 2>&1 \
  || { echo "refused: the screen's run commit $SCREENED is not in the history of $COMMIT" >&2; finish 84; }
git diff --quiet "$SCREENED" "$COMMIT" -- $SCREEN_FILES \
  || { echo "refused: a file the screen read differs between its run commit $SCREENED and $COMMIT" >&2; finish 85; }
O=runs/hearth_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$FORMS" > "$S/forms"
echo "$SCREENED" > "$S/screened"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_trend_line_h.py tests/test_heart_trend_line_h_launcher.py tests/test_h35_heart_trend_forms.py \
  tests/test_run_h35_launcher.py tests/test_heart_trend_line_g.py tests/test_heart_trend_line_g_launcher.py \
  tests/test_heart_trend_line_f.py tests/test_heart_trend_line_d.py tests/test_heart_trend_line.py \
  tests/test_build_heart_organ_scope.py tests/test_heart_trend_build.py tests/test_b0_trend.py \
  tests/test_h28_heart_trend_strength.py tests/test_halve.py tests/test_hidden_counts.py tests/test_s2_rules.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, its half, `vec check` on both; exit 11 build, 12 halve, 13 check,
           # 14 no stamp, 15 stamp not this commit and clean
  local d="$O/$1" name=$2 st
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py heart "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  $PY experiments/builds/halve.py "$d/$name.h5ad" "$d/${name}_half.h5ad" > /dev/null 2>&1 || return 12
  for f in "$d/$name.h5ad" "$d/${name}_half.h5ad"; do
    uv run --frozen --no-sync vec check "$f" --board "$BOARD" > /dev/null 2>&1 || return 13
  done
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
identity() {  # step name args...: the step's build, then both files' sha256 against the recorded ones
  local name=$1 want=$2 want_half=$3 code=$4 code_half=$5 got
  shift 5
  echo running > "$S/$STEP.exit"
  build "$STEP" "$name" "$@" || finish $?
  got=$(sha256sum "$O/$STEP/$name.h5ad" | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "stopped: $STEP's command does not give its recorded file" >&2; finish "$code"; }
  got=$(sha256sum "$O/$STEP/${name}_half.h5ad" | cut -d' ' -f1)
  [ "$got" = "$want_half" ] || { echo "stopped: $STEP's half is not its recorded half" >&2; finish "$code_half"; }
  step_ok
}
: > "$S/stamps"

STEP=D
identity T2heart_organ "$D_SHA" "$D_HALF_SHA" 20 21 $D_ARGS
STEP=OTS05
identity T2heart_organ_ts05 "$OTS05_SHA" "$OTS05_HALF_SHA" 22 23 $OTS05_ARGS
STEP=OTS10
identity T2heart_organ_ts10 "$OTS10_SHA" "$OTS10_HALF_SHA" 24 25 $OTS10_ARGS
STEP=OTS075
identity T2heart_organ_ts075 "$OTS075_SHA" "$OTS075_HALF_SHA" 26 27 $OTS075_ARGS
STEP=OTS125
identity T2heart_organ_ts125 "$OTS125_SHA" "$OTS125_HALF_SHA" 28 29 $OTS125_ARGS
STEP=OTS15
identity T2heart_organ_ts15 "$OTS15_SHA" "$OTS15_HALF_SHA" 32 33 $OTS15_ARGS

TAGS=""
for form in $FORMS; do
  case "$form" in
    category) tag=tk; args=$TK_ARGS ;;
    contrast) tag=tx; args=$TX_ARGS ;;
    whole) tag=tw; args=$TW_ARGS ;;
    portion) tag=tp; args=$TP_ARGS ;;
  esac
  TAGS="$TAGS $tag"
  for STEP in "${tag}a" "${tag}b"; do
    echo running > "$S/$STEP.exit"
    build "$STEP" "T2heart_organ_ts125_$tag" $args || finish $?
    step_ok
  done
done

STEP=done
echo running > "$S/$STEP.exit"
dirs="D OTS05 OTS10 OTS075 OTS125 OTS15"
for tag in $TAGS; do dirs="$dirs ${tag}a ${tag}b"; done
(cd "$O" && for d in $dirs; do sha256sum "$d"/*.h5ad || exit 1; done) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for tag in $TAGS; do
  for f in "T2heart_organ_ts125_$tag" "T2heart_organ_ts125_${tag}_half"; do
    cmp -s "$O/${tag}a/$f.h5ad" "$O/${tag}b/$f.h5ad" || finish 31
  done
done
step_ok
