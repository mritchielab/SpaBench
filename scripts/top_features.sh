#!/usr/bin/env bash
# usage: ./top_features.sh /path/to/root
ROOT="${1:-.}"
SIZES=(500 1000 2000 5000 10000)

find "$ROOT" -type f -name '04_*.csv' -print0 | while IFS= read -r -d '' f; do
  base=$(basename "$f")
  # skip files like 04_xxx_2000.csv
  if [[ "$base" =~ _[0-9]+\.csv$ ]]; then
    continue
  fi
  dirn=$(dirname "$f")
  echo "Processing: $f"
  lines=$(wc -l < "$f" | tr -d '[:space:]')
  for n in "${SIZES[@]}"; do
    out="$dirn/top_${n}.csv"
    [ -e "$out" ] && { echo "  Exists, skip: $out"; continue; }
    if [ "$lines" -ge "$n" ]; then
      head -n "$n" "$f" > "$out" && echo "  Wrote top $n -> $out"
    else
      cp "$f" "$out" && echo "  Copied full (lines=$lines) -> $out" || echo "  ERROR copying -> $out"
    fi
  done
done

