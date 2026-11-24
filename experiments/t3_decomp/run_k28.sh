#!/usr/bin/env bash
# Row K28 (experiments/harness/TRIES.md; experiments/t3_decomp/k28_chip_partners.py, T3 candidate CHIPQ-C) on a DGX
# campaign checkout (never the reference checkout), the machine that holds the pinned database. One step per call. Run
# from the checkout's root, detached (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a
# foreground sleep.
#
#   bash experiments/t3_decomp/run_k28.sh g0a COMMIT      G0a, the first step of the database-only gate G0
#
# COMMIT: the full hash of the commit the step runs at, the row's registration commit (row K28 in TRIES.md, the
# amendment of register row X018 in EXTERNAL_DATA.md, G0a's code and tests). The steps after it (G0b, G0c, the read) are
# added here with their code, which by the specification is written only after G0a has passed on the real database.
#
# Fixed folders (no override): the fetch ~/v/vec-ext/t3/chipatlas (chipq_fetch.DATA_ROOT, pinned on 2026-11-23); the
# scratch folder S = ~/v/scratch/k28/<first 7 of COMMIT>/g0a (TMPDIR under it). G0 reads database files only: no data
# folder, no run folder and no table folder is named here.
#
# g0a, stopping at the first failure:
#   1. HEAD is COMMIT (exit 97), the tracked tree is clean (98), row K28 is in TRIES.md (95); chipq_kept.yaml and
#      CHIPQ_KEPT.sha256 are tracked, the fetch folder is a real directory (no symlink), and EXTERNAL_DATA.md's row
#      X018 carries the amendment's opening words (96);
#   2. a reading is in the row's results folder already (the candidate is closed): made at this HEAD, 28 again;
#      made at another commit, naming none, or unreadable, 44. Nothing runs;
#   3. the targeted tests, hand-made files only, 8 BLAS threads (10);
#   4. G0a (k28_chip_partners.py g0a --commit COMMIT; it hashes every file of the fetch manifest again before it
#      counts): 0 passes; 3 a bar is not met: the step has made the STOP reading, once (28: no fold read, no try, no
#      file; 41 if the reading is not there); 2 refused (22); other 20. Its three lines of counts go to $S/g0a.out and
#      to this launcher's own output.
#   all.exit: 0, 10, 20, 22, 28, 41, 44, 95, 96, 97, 98; 99 the launcher stopped unexpectedly.
#
# Every step writes its exit code to $S/<step>.exit and every path writes $S/all.exit. A call without the right
# arguments writes $HOME/v/scratch/k28/usage.exit (2). Every process's output goes to /dev/null, except the gate's
# three lines of counts. Nothing is ever deleted, fetched or uploaded. NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under $S; 8
# BLAS threads; PATH starts with $HOME/.local/bin before the first uv call.
set -u
STEP=${1:-}
if [ "$STEP" = g0a ] && [ "$#" -eq 2 ] && [[ "$2" =~ ^[0-9a-f]{40}$ ]]; then
  COMMIT=$2
else
  mkdir -p "$HOME/v/scratch/k28" && echo 2 > "$HOME/v/scratch/k28/usage.exit"
  echo "usage: run_k28.sh g0a COMMIT (the full 40-hex hash of the commit the step runs at)" >&2
  exit 2
fi
S="$HOME/v/scratch/k28/${COMMIT:0:7}/$STEP"
mkdir -p "$S/tmp" || exit 99
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/all.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT

export PATH="$HOME/.local/bin:$PATH"
CHIPATLAS=$HOME/v/vec-ext/t3/chipatlas  # the fetch (chipq_fetch.DATA_ROOT); fixed, no override
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
SCRIPT=experiments/t3_decomp/k28_chip_partners.py
OUT=experiments/t3_decomp/results/k28_chip_partners
WORDS="Reopened by the captain's order of 2026-11-24"  # the opening words of row X018's amendment
HEAD_NOW=$(git rev-parse HEAD 2> /dev/null) || HEAD_NOW=""

# ------------------------------------------------------------------------------------------------------ the gates
[ "$HEAD_NOW" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 97; }
[ -z "$(git status --porcelain --untracked-files=no 2>&1)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
grep -q '^| K28 ' experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: no TRIES.md row K28" >&2; finish 95; }
for f in experiments/t3_ext/chipq_kept.yaml experiments/t3_ext/CHIPQ_KEPT.sha256; do
  GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$f" > /dev/null 2>&1 \
    || { echo "refused: $f is not tracked" >&2; finish 96; }
done
{ [ -d "$CHIPATLAS" ] && [ ! -L "$CHIPATLAS" ]; } \
  || { echo "refused: $CHIPATLAS is not a real directory" >&2; finish 96; }
grep '^| X018 ' EXTERNAL_DATA.md 2> /dev/null | grep -qF "$WORDS" \
  || { echo "refused: row X018 of EXTERNAL_DATA.md does not carry the opening words of its amendment" >&2; finish 96; }
git rev-parse HEAD > "$S/head"

# a reading is there already: the candidate is closed, and nothing runs again
if [ -e "$OUT/reading.json" ] || [ -L "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ] \
    || [ -L "$OUT/reading.txt" ]; then
  got=$($PY -c 'import json, sys; r = json.load(open(sys.argv[1])); print("same" if (r.get("reading_code") or {}).get("commit") == sys.argv[2] else "other")' \
    "$OUT/reading.json" "$HEAD_NOW" 2> /dev/null) || got=unreadable
  if [ "$got" = same ]; then finish 28; fi
  finish 44
fi

# ------------------------------------------------------------------------------------------------------ the tests
threads 8
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_chipq_partners.py tests/test_k28_chip_partners.py tests/test_k28_launcher.py tests/test_chipq.py \
  tests/test_k27_chip_query.py > /dev/null 2>&1 || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

# ------------------------------------------------------------------------------------------------------ G0a
rc=0
$PY "$SCRIPT" g0a --commit "$COMMIT" > "$S/g0a.out" 2> /dev/null || rc=$?
echo $rc > "$S/g0a.exit"
if [ $rc = 0 ] || [ $rc = 3 ]; then cat "$S/g0a.out"; fi
if [ $rc = 3 ]; then
  if [ -e "$OUT/reading.json" ] && [ -e "$OUT/reading.txt" ]; then finish 28; fi
  finish 41
fi
if [ $rc = 2 ]; then finish 22; fi
if [ $rc != 0 ]; then finish 20; fi
finish 0
