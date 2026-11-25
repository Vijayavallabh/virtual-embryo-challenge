#!/usr/bin/env bash
# Parallel range download of the MOSTA files pinned in experiments/heart_ext/fetch_oic.py (register row X022), for a
# server that caps each connection (about 0.2 MB/s each on 2026-11-25): each pinned file is split into 32 MiB ranges,
# fetched in parallel with curl, each range checked for its exact length, then joined in order into <out>/<name>,
# whose length must equal the pin and which must pass `gzip -t` (the gzip CRC-32 of the whole file). fetch_oic.py,
# run after this, finds each complete file of its pinned size, records its sha256 in MANIFEST.json and fetches the
# Xenium member itself. Only the pinned URLs are requested (read from fetch_oic.py's PINNED); nothing is deleted.
#   bash fetch_oic_chunks.sh <path of fetch_oic.py> <out folder> <chunk folder> <parallel connections>
set -euo pipefail
PY=${1:?usage: fetch_oic_chunks.sh FETCH_OIC_PY OUT CHUNKS PARALLEL}
OUT=${2:?usage: fetch_oic_chunks.sh FETCH_OIC_PY OUT CHUNKS PARALLEL}
CH=${3:?usage: fetch_oic_chunks.sh FETCH_OIC_PY OUT CHUNKS PARALLEL}
P=${4:?usage: fetch_oic_chunks.sh FETCH_OIC_PY OUT CHUNKS PARALLEL}
CHUNK=$((32 * 1024 * 1024))
mkdir -p "$OUT" "$CH"
python3 -I -B -c '
import importlib.util, sys
spec = importlib.util.spec_from_file_location("fetch_oic_pinned", sys.argv[1])
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
for e in m.PINNED:
    if e["kind"] == "file":
        assert e["url"].startswith("https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/Bin1_matrix/"), e["url"]
        print(e["name"], e["url"], e["bytes"], sep="\t")
' "$PY" > "$CH/list.tsv"
[ "$(wc -l < "$CH/list.tsv")" = 23 ] || { echo "refused: the pinned list does not hold 23 MOSTA files" >&2; exit 2; }

fetch_chunk() {  # URL START END FILE: one range, retried until it has exactly END - START + 1 bytes
  local url=$1 s=$2 e=$3 f=$4 want n=0
  want=$((e - s + 1))
  while :; do
    if [ -f "$f" ] && [ "$(stat -c %s "$f")" = "$want" ]; then return 0; fi
    n=$((n + 1))
    [ "$n" -le 80 ] || { echo "FAILED $f" >&2; return 1; }
    if curl -sS -f -A "curl (public data download)" --connect-timeout 30 --max-time 1800 -r "$s-$e" \
         -o "$f.tmp" "$url"; then
      mv -f "$f.tmp" "$f"
    else
      sleep 5
    fi
  done
}
export -f fetch_chunk

jobs_file="$CH/jobs.txt"
: > "$jobs_file"
while IFS=$'\t' read -r name url size; do
  [ -f "$OUT/$name" ] && continue
  mkdir -p "$CH/$name"
  n=$(( (size + CHUNK - 1) / CHUNK ))
  for ((i = 0; i < n; i++)); do
    s=$((i * CHUNK)); e=$((s + CHUNK - 1)); [ "$e" -lt "$size" ] || e=$((size - 1))
    printf '%s %s %s %s\n' "$url" "$s" "$e" "$CH/$name/$(printf '%06d' "$i")" >> "$jobs_file"
  done
done < "$CH/list.tsv"
echo "chunks to fetch: $(wc -l < "$jobs_file")"
xargs -P "$P" -n 4 bash -c 'fetch_chunk "$@"' _ < "$jobs_file"

while IFS=$'\t' read -r name url size; do
  [ -f "$OUT/$name" ] && continue
  n=$(( (size + CHUNK - 1) / CHUNK ))
  : > "$OUT/$name.joining"
  for ((i = 0; i < n; i++)); do cat "$CH/$name/$(printf '%06d' "$i")" >> "$OUT/$name.joining"; done
  [ "$(stat -c %s "$OUT/$name.joining")" = "$size" ] || { echo "refused: $name joined to the wrong length" >&2; exit 2; }
  gzip -t "$OUT/$name.joining" || { echo "refused: $name fails gzip -t" >&2; exit 2; }
  mv "$OUT/$name.joining" "$OUT/$name"
  echo "joined $name $size"
done < "$CH/list.tsv"
echo "all 23 MOSTA files joined"
