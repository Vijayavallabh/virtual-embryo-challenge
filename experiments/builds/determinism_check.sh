#!/bin/bash
# Determinism check of the ten test-board builds (docs/PLAN.md decision log 2026-11-05 (b)). It rebuilds every set
# (s1, s2, t3alt) with experiments/builds/build_all.sh at each BLAS thread count given, each into its own folder, then
# compares every file's SHA256 with an archived set. It protects the freeze-morning hash check: a file that changes
# with the thread count, or between runs, would make that check fail for a reason that is not a method change.
# Target-free: the same builders on the released inputs; it reads no target and computes no score. Nothing is
# uploaded.
#
#   bash experiments/builds/determinism_check.sh ARCHIVE_SHA256SUMS THREADS...
#   e.g. bash experiments/builds/determinism_check.sh runs/archive/builds_ae6d6b3/SHA256SUMS 1 4 8 16
#
# Run on the DGX from the repo root at a clean, committed tree, with no sync while it runs. The thread counts run in
# parallel, one chain each; the outputs go to runs/determinism_<HEAD>/t<threads>/ (ignored by git), logs next to them.
set -euo pipefail
sums=${1:?usage: determinism_check.sh ARCHIVE_SHA256SUMS THREADS...}
shift
if [ $# -lt 1 ]; then
  echo "give at least one thread count" >&2
  exit 2
fi
sums=$(readlink -f "$sums")
head=$(git rev-parse --short HEAD)
root=runs/determinism_$head
if [ -e "$root" ]; then
  echo "refused: $root exists" >&2
  exit 2
fi
mkdir -p "$root"
pids=()
for t in "$@"; do
  (
    for set_name in s1 s2 t3alt; do
      BUILD_OUT="$root/t$t" BUILD_THREADS="$t" bash experiments/builds/build_all.sh "$set_name"
    done
  ) >"$root/t$t.log" 2>&1 &
  pids+=("$!")
done
failed=0
for p in "${pids[@]}"; do
  wait "$p" || failed=1
done
if [ "$failed" -ne 0 ]; then
  echo "a build chain failed; see $root/t*.log" >&2
  exit 1
fi
status=0
for t in "$@"; do
  echo "== BLAS threads $t, against $sums"
  (cd "$root/t$t" && sha256sum -c "$sums") || status=1
done
if [ "$status" -eq 0 ]; then
  echo "all files byte-identical to the archived set at every thread count"
else
  echo "MISMATCH: at least one file differs from the archived set" >&2
fi
exit "$status"
