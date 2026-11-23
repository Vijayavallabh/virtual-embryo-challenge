#!/usr/bin/env bash
# Row K27 (experiments/harness/TRIES.md; experiments/t3_decomp/k27_chip_query.py, T3 candidate CHIPQ-S) on a DGX campaign
# checkout (~/v/vec-campaign-c, never the reference checkout), the machine rows K18 and K21 to K24 ran on, in TWO
# invocations. Run from the checkout's root, detached (setsid nohup ... &), and wait on $S/all.exit with a background
# poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/run_k27.sh db REG                 the database part, at the row's registration commit
#   bash experiments/t3_decomp/run_k27.sh run SECOND REG         the counted part, at the row's second commit
#   Before the database part runs, the folder ~/v/vec-ext/t3/tables_k27 must exist (mkdir -p, outside every checkout).
#
# REG: the full hash of the row's registration commit (the row is in TRIES.md; TABLES_K27.sha256 is tracked and empty;
# the 11 table pins and K27_LISTING_SHA256 in the script are None). SECOND: the full hash of the row's second commit,
# made after P2 and before P3, which fills the listing's 11 lines and the pins; REG is its ancestor, and REG..SECOND
# changes only experiments/t3_ext/TABLES_K27.sha256 and experiments/t3_decomp/k27_chip_query.py. Fixed folders (no
# override): the K27 tables ~/v/vec-ext/t3/tables_k27 (VEC_T3_EXT_TABLES_K27), `two` ~/v/vec-ext/t3/tables_two_study
# (VEC_T3_EXT_TABLES_TWO), the fetch ~/v/vec-ext/t3/chipatlas; the data $PWD/data/raw; the reference checkout's gata4 run
# folder. The run folder R = runs/k27_<first 7 of REG>, the same for both invocations (resumable); the scratch folder
# S = ~/v/scratch/k27_<first 7 of REG>/db or /run (TMPDIR under it). K23's run folder runs/k23_<first 7 of K23's commit>
# must be in this checkout (guard 2 reads its per-seed x_e2 at the sha256 K23's reading records).
#
# Invocation 1, db (nothing counted), stopping at the first failure:
#   1. HEAD is REG (exit 97), the tracked tree is clean (98), the row naming the script is in TRIES.md (95);
#      TABLES_K27.sha256 is tracked and empty and the script's 11 table pins and its listing pin are None (22 otherwise);
#      chipq_kept.yaml, CHIPQ_KEPT.sha256 and K23's committed pregate.json and reading.json are tracked, the K27 tables
#      folder and the fetch folder are real directories (no symlink), EXTERNAL_DATA.md holds row X018 naming ChIP-Atlas
#      (96);
#   2. the targeted tests, synthetic data only, 8 BLAS threads (10);
#   3. chipq_fetch.py verify: every manifest file re-hashed (its exit 2, or any other non-zero: 22);
#   4. P1 (chipq.py p1, Ctnnb1's two integers into R/p1.out, its exit into R/p1.exit, nothing else written): 0 passes,
#      3 fails (the STOP reading is made once, then 28: no try, no file, no second commit), 2 refused (22), other 20;
#   5. P2 (k27_chip_query.py p2 --out R; skipped when R/p2.json exists, whose "passes" is then read): 0 passes (the 11
#      tables written into the K27 folder), 3 fails (the STOP reading once, then 29), 2 refused (22: pins filled, a
#      placebo table equal to the query's, a sha256 twice, a scaled SD off S...), other 20;
#   6. the second commit's lines (k27_chip_query.py pins --runs R) into $S/second_commit.txt: the 11 listing lines,
#      K27_LISTING_SHA256 and K27_TABLE_SHA256 (sha256 and table names only; its refusal 22).
#   all.exit: 0 (P2 passed, the lines written), 10, 20, 22, 28, 29, 41, 43, 44 (the STOP reading crashed, was
#   INCOMPLETE, or found another commit's reading), 95, 96, 97, 98; 99 the launcher stopped unexpectedly.
#
# Invocation 2, run (the eight tries count from the start of step 5, P3):
#   1. HEAD is SECOND and REG is a strict ancestor of it (97), the tracked tree is clean (98), the row is in TRIES.md
#      (95); the listing has its 11 lines and the 11 pins and the listing pin are filled, and REG..SECOND changes only
#      the listing and the script (22); invocation 1 finished 0 at REG ($HOME/v/scratch/k27_<REG7>/db/all.exit) and
#      R/p2.json is there, the tracked inputs as in invocation 1, all three outside folders real directories, X018 (96);
#   2. the targeted tests (10);
#   3. chipq_fetch.py verify (22);
#   4. the tables checked against the listing, the pins, their manifests and the gene order (k27_chip_query.py check:
#      2 -> 22, other -> 20; nothing counted);
#   5. P3 (k27_chip_query.py pregate --out R; skipped when R/pregate.json exists): 0, 5 the arm is dropped on either
#      correlation (the drop's reading once, then 25; the eight tries counted), 2 refused (22), other 20; a rerun whose
#      pregate.json records a drop finishes 25 again without a second reading;
#   6. the GUARD PART: barred_form and x_e2 only, both worlds, ten seeds, five parallel parts into R/guard/ (30 if a part
#      did not exit 0), then the guard check (R/guard.json): its 6 makes the VOID reading once and finishes 26; 2 -> 22,
#      3 -> 30, other 20; a rerun over a guard.json that does not pass finishes 26 again with no scoring part started;
#   7. the five scoring parts in parallel: plain g1 and g2, each as seeds 0-4 and 5-9 (8 BLAS threads a process), and
#      the depth-matched world (8 workers by 2 threads, 16 BLAS threads in its set-up); each part's exit code in
#      R/exit/<part>, its log in R/logs/<part>.log (never printed); log lines other than the expected kinds are counted,
#      never printed (30 if a part did not exit 0);
#   8. the read into experiments/t3_decomp/results/k27_chip_query/ (output to /dev/null): 0 a reading (VOID included),
#      43 INCOMPLETE, 44 a reading is there already (never overwritten), 41 a crash or a refusal in the read.
#   all.exit: 0, 10, 20, 22, 25, 26, 30, 41, 43, 44, 95, 96, 97, 98; 99 the launcher stopped unexpectedly.
#
# A stop's reading (28, 29, 25, 26) is made once: when OUT/reading.json exists and its reading_code.commit is this HEAD
# the step finishes with its code without reading again; one made at another commit, naming none, or unreadable
# finishes 44. Every step writes its exit code to $S/<step>.exit and every path writes $S/all.exit. A call without the
# right arguments writes $HOME/v/scratch/k27/usage.exit (2). Nothing is ever removed. NUMPY_MADVISE_HUGEPAGE=0; TMPDIR
# under $S; at most 16 BLAS threads a process; PATH starts with $HOME/.local/bin before the first uv call.
set -u
MODE=${1:-}
if [ "$MODE" = db ] && [ "$#" -eq 2 ] && [ -n "$2" ]; then
  REG=$2
  SECOND=""
elif [ "$MODE" = run ] && [ "$#" -eq 3 ] && [ -n "$2" ] && [ -n "$3" ]; then
  SECOND=$2
  REG=$3
else
  mkdir -p "$HOME/v/scratch/k27" && echo 2 > "$HOME/v/scratch/k27/usage.exit"
  echo "usage: run_k27.sh db REGISTRATION_COMMIT | run_k27.sh run SECOND_COMMIT REGISTRATION_COMMIT" >&2
  exit 2
fi
S="$HOME/v/scratch/k27_${REG:0:7}/$MODE"
mkdir -p "$S/tmp" || exit 99
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/all.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT

export PATH="$HOME/.local/bin:$PATH"
TABLES_K27=$HOME/v/vec-ext/t3/tables_k27  # the K27 tables (P2 writes them); fixed, no override
TABLES_TWO=$HOME/v/vec-ext/t3/tables_two_study  # `two`, K22's and K23's; fixed, no override
CHIPATLAS=$HOME/v/vec-ext/t3/chipatlas  # the fetch (chipq_fetch.DATA_ROOT); fixed, no override
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv" \
  VEC_T3_EXT_TABLES_K27="$TABLES_K27" VEC_T3_EXT_TABLES_TWO="$TABLES_TWO"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
SCRIPT=experiments/t3_decomp/k27_chip_query.py
CHIPQ=experiments/t3_ext/chipq.py
FETCH=experiments/t3_ext/chipq_fetch.py
LISTING=experiments/t3_ext/TABLES_K27.sha256
DATA="$PWD/data/raw"  # this checkout's own data, never the installed package's
R=runs/k27_${REG:0:7}
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
OUT=experiments/t3_decomp/results/k27_chip_query
HEAD_NOW=$(git rev-parse HEAD 2> /dev/null) || HEAD_NOW=""

# ------------------------------------------------------------------------------------------------------ the gates
if [ "$MODE" = db ]; then
  [ "$HEAD_NOW" = "$REG" ] || { echo "refused: HEAD is not $REG" >&2; finish 97; }
else
  [ "$HEAD_NOW" = "$SECOND" ] || { echo "refused: HEAD is not $SECOND" >&2; finish 97; }
  { [ "$REG" != "$SECOND" ] && git merge-base --is-ancestor "$REG" "$SECOND" 2> /dev/null; } \
    || { echo "refused: $REG is not a strict ancestor of $SECOND" >&2; finish 97; }
fi
[ -z "$(git status --porcelain --untracked-files=no 2>&1)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
grep -q '^| .*`experiments/t3_decomp/k27_chip_query.py`' experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: no TRIES.md row names the script" >&2; finish 95; }
GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$LISTING" > /dev/null 2>&1 \
  || { echo "refused: $LISTING is not tracked" >&2; finish 22; }
# the pin lines: the listing pin and the 11 table pins, exactly as the script writes them
n_none=$(grep -cE '^    "chipq_(Gata4|rand_seed[0-9])": None,$' "$SCRIPT" 2> /dev/null) || n_none=0
n_set=$(grep -cE '^    "chipq_(Gata4|rand_seed[0-9])": "[0-9a-f]{64}",$' "$SCRIPT" 2> /dev/null) || n_set=0
l_none=$(grep -cE '^K27_LISTING_SHA256: str \| None = None$' "$SCRIPT" 2> /dev/null) || l_none=0
l_set=$(grep -cE '^K27_LISTING_SHA256: str \| None = "[0-9a-f]{64}"$' "$SCRIPT" 2> /dev/null) || l_set=0
n_lines=$(grep -cE '^[0-9a-f]{64}  chipq_(Gata4|rand_seed[0-9])$' "$LISTING" 2> /dev/null) || n_lines=0
if [ "$MODE" = db ]; then
  { [ "$n_none" = 11 ] && [ "$n_set" = 0 ] && [ "$l_none" = 1 ] && [ "$l_set" = 0 ] && [ -f "$LISTING" ] \
      && [ ! -s "$LISTING" ]; } \
    || { echo "refused: the pins or the listing are not empty (the db part runs at the registration commit)" >&2;
         finish 22; }
else
  { [ "$n_set" = 11 ] && [ "$n_none" = 0 ] && [ "$l_set" = 1 ] && [ "$l_none" = 0 ] && [ "$n_lines" = 11 ] \
      && [ "$(wc -l < "$LISTING")" = 11 ]; } \
    || { echo "refused: the pins or the listing are not filled (the run part runs at the second commit)" >&2;
         finish 22; }
  changed=$(git diff --name-only "$REG" "$SECOND" 2> /dev/null | sort | tr '\n' ' ') || changed=unreadable
  [ "$changed" = "experiments/t3_decomp/k27_chip_query.py experiments/t3_ext/TABLES_K27.sha256 " ] \
    || { echo "refused: the second commit changes more than the listing and the pins" >&2; finish 22; }
fi
for f in experiments/t3_ext/chipq_kept.yaml experiments/t3_ext/CHIPQ_KEPT.sha256 \
    experiments/t3_decomp/results/k23_outside_exponents/pregate.json \
    experiments/t3_decomp/results/k23_outside_exponents/reading.json; do
  GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$f" > /dev/null 2>&1 \
    || { echo "refused: $f is not tracked" >&2; finish 96; }
done
DIRS=("$TABLES_K27" "$CHIPATLAS")
if [ "$MODE" = run ]; then DIRS+=("$TABLES_TWO"); fi
for d in "${DIRS[@]}"; do
  { [ -d "$d" ] && [ ! -L "$d" ]; } || { echo "refused: $d is not a real directory" >&2; finish 96; }
done
grep '^| X018 ' EXTERNAL_DATA.md 2> /dev/null | grep -qF ChIP-Atlas \
  || { echo "refused: EXTERNAL_DATA.md holds no row X018 naming ChIP-Atlas" >&2; finish 96; }
if [ "$MODE" = run ]; then
  [ "$(cat "$HOME/v/scratch/k27_${REG:0:7}/db/all.exit" 2> /dev/null)" = 0 ] \
    || { echo "refused: the database part has not finished 0 at $REG" >&2; finish 96; }
  [ -f "$R/p2.json" ] || { echo "refused: $R/p2.json is missing" >&2; finish 96; }
fi
git rev-parse HEAD > "$S/head"
mkdir -p "$R/exit" "$R/logs"

# ------------------------------------------------------------------------------------------------------ helpers
# a stop's reading (28, 29, 25, 26), made once: this HEAD's reading there already -> the same finish; another -> 44
stop_read() {
  local rc=0 got
  if [ -e "$OUT/reading.json" ]; then
    got=$($PY -c 'import json, sys; r = json.load(open(sys.argv[1])); print("same" if (r.get("reading_code") or {}).get("commit") == sys.argv[2] else "other")' \
      "$OUT/reading.json" "$HEAD_NOW" 2> /dev/null) || got=unreadable
    if [ "$got" = same ]; then finish "$1"; fi
    finish 44
  fi
  $PY "$SCRIPT" read --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
  echo $rc > "$S/stop_read.exit"
  if [ $rc = 3 ]; then finish 43; fi
  if [ $rc = 4 ]; then finish 44; fi
  if [ $rc != 0 ]; then finish 41; fi
  finish "$1"
}

json_key() {  # $1 = a json file, $2 = a python expression on r: prints its value, or "unreadable"
  $PY -c 'import json, sys; r = json.load(open(sys.argv[1])); print(eval(sys.argv[2], {"r": r}))' "$1" "$2" \
    2> /dev/null || echo unreadable
}

threads 8
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_k27_chip_query.py tests/test_k27_launcher.py tests/test_chipq.py tests/test_k24_wider_tables.py \
  tests/test_build_k24_form.py tests/test_k23_outside_exponents.py tests/test_k22_outside_forms.py \
  tests/test_t3_outside_forms.py tests/test_k21_factor_forms.py tests/test_k18_built_files.py > "$S/tests.log" 2>&1 \
  || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

rc=0
$PY "$FETCH" verify > "$S/verify.log" 2>&1 || rc=$?
echo $rc > "$S/verify.exit"
if [ $rc != 0 ]; then finish 22; fi

# ------------------------------------------------------------------------------------------------------ invocation 1
if [ "$MODE" = db ]; then
  rc=0
  $PY "$CHIPQ" p1 --data "$DATA" --run-dir "$RUN_DIR" > "$R/p1.out" 2> "$S/p1.log" || rc=$?
  echo $rc > "$S/p1.exit"
  echo $rc > "$R/p1.exit"
  if [ $rc = 3 ]; then stop_read 28; fi
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc != 0 ]; then finish 20; fi
  if [ ! -e "$R/p2.json" ]; then
    rc=0
    $PY "$SCRIPT" p2 --out "$R" --data "$DATA" --run-dir "$RUN_DIR" > "$S/p2.log" 2>&1 || rc=$?
    echo $rc > "$S/p2.exit"
    if [ $rc = 3 ]; then stop_read 29; fi
    if [ $rc = 2 ]; then finish 22; fi
    if [ $rc != 0 ]; then finish 20; fi
  fi
  p2pass=$(json_key "$R/p2.json" 'r.get("passes") is True')
  if [ "$p2pass" = False ]; then stop_read 29; fi
  if [ "$p2pass" != True ]; then finish 20; fi
  rc=0
  $PY "$SCRIPT" pins --runs "$R" > "$S/second_commit.txt" 2> "$S/pins.log" || rc=$?
  echo $rc > "$S/pins.exit"
  if [ $rc != 0 ]; then finish 22; fi
  finish 0
fi

# ------------------------------------------------------------------------------------------------------ invocation 2
rc=0
$PY "$SCRIPT" check --runs "$R" --data "$DATA" --run-dir "$RUN_DIR" > "$S/check.log" 2>&1 || rc=$?
echo $rc > "$S/check.exit"
if [ $rc = 2 ]; then finish 22; fi
if [ $rc != 0 ]; then finish 20; fi

if [ ! -e "$R/pregate.json" ]; then
  rc=0
  $PY "$SCRIPT" pregate --out "$R" --data "$DATA" --run-dir "$RUN_DIR" > "$S/pregate.log" 2>&1 || rc=$?
  echo $rc > "$S/pregate.exit"
  if [ $rc = 5 ]; then stop_read 25; fi
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc != 0 ]; then finish 20; fi
fi
pdrop=$(json_key "$R/pregate.json" '(r.get("pregate") or {}).get("dropped") is not None')
if [ "$pdrop" = True ]; then stop_read 25; fi
if [ "$pdrop" != False ]; then finish 20; fi

OK_LINE='^(K27( guard)? plain g[12] seed [0-9]+: (scored|present, skipped)|wrote [^ ]+\.json.*|\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s|self-check passed .*)$'
PARTS=("plain_g1_0-4" "plain_g2_0-4" "plain_g1_5-9" "plain_g2_5-9" "dm")

parts() {  # $1 = guard|score: the five parts at once; returns 1 when a part did not exit 0
  local pids=() part p bad=0 n pre
  pre=""
  if [ "$1" = guard ]; then pre="guard_"; fi
  for part in "${PARTS[@]}"; do
    (
      if [ "$part" = dm ]; then
        threads 16
        $PY "$SCRIPT" run --part "$1" --world dm --workers 8 --threads 2 --out "$R" --data "$DATA" \
          --run-dir "$RUN_DIR" > "$R/logs/$pre$part.log" 2>&1
      else
        threads 8
        IFS=_ read -r _ task seeds <<< "$part"
        $PY "$SCRIPT" run --part "$1" --world plain --task "$task" --seeds "$seeds" --out "$R" --data "$DATA" \
          --run-dir "$RUN_DIR" > "$R/logs/$pre$part.log" 2>&1
      fi
      echo $? > "$R/exit/$pre$part"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || true; done
  for part in "${PARTS[@]}"; do
    [ "$(cat "$R/exit/$pre$part" 2> /dev/null)" = 0 ] || bad=1
    n=$(grep -cvE "$OK_LINE" "$R/logs/$pre$part.log" 2> /dev/null) || true
    echo "${n:-unreadable}" > "$R/exit/$pre$part.other_lines"
  done
  return $bad
}

if [ ! -e "$R/guard.json" ]; then
  threads 8
  parts guard || finish 30
  rc=0
  $PY "$SCRIPT" guard --runs "$R" > "$S/guard.log" 2>&1 || rc=$?
  echo $rc > "$S/guard.exit"
  if [ $rc = 6 ]; then stop_read 26; fi
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc = 3 ]; then finish 30; fi
  if [ $rc != 0 ]; then finish 20; fi
fi
# a guard.json that does not pass (a rerun after a guard failure) finishes 26 again; no scoring part starts
gpass=$(json_key "$R/guard.json" 'r.get("passes") is True')
if [ "$gpass" = False ]; then stop_read 26; fi
if [ "$gpass" != True ]; then finish 20; fi

parts score || finish 30

rc=0
$PY "$SCRIPT" read --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
echo $rc > "$S/read.exit"
if [ $rc = 3 ]; then finish 43; fi
if [ $rc = 4 ]; then finish 44; fi
if [ $rc != 0 ]; then finish 41; fi
finish 0
